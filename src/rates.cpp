#include "orpp/csv.hpp"
#include <limits.h>
#include "districts.hpp"

using namespace std;
using namespace orpp;

constexpr int maxint = numeric_limits<int>::max();

enum cohorts { c0, c20, c65, c80, numcohorts};

cohorts v2cohort(unsigned v)
{
    if(v < 20)
        return c0;
    else
    {
        if(v<65)
            return c20;
        else
        {
            if(v < 80)
                return c65;
            else
            {
                if(v < 105)
                    return c80;
                else
                    return numcohorts;
            }
        }
    }
};

struct counter
{
    unsigned none = 0;
    unsigned wrong = 0;
    unsigned under = 0;
    unsigned over = 0;
};

int date2int(const string s)
{
    struct tm ti;
    ti.tm_sec=0;
    ti.tm_min=0;
    ti.tm_hour=0;
    ti.tm_mday= stoi(s.substr(8,2));
    ti.tm_mon = stoi(s.substr(5,2))-1;
    ti.tm_year=stoi(s.substr(0,4))-1900;

    time_t res = mktime(&ti);
    int di = round((res + 86400 *0.5) / 86400);
    return di;
}

int date2int(const string s, int zerodate, int lastdate, counter& c)
{
    if(s=="")
    {
        c.none++;
        return maxint;
    }
    if(s.size() != 10)
    {
        c.wrong++;
        return maxint;
    }
    int di = date2int(s);
    if(di < zerodate)
    {
        c.under++;
        return maxint;
    }
    else if(di > lastdate)
    {
        c.over++;
        return maxint;
    }
    else
        return di;
}

int str2district(const string& okres_lau_kod)
{
   for(unsigned i=0; i<districts.size(); i++)
   {
       if(okres_lau_kod==districts[i].code)
           return i;
   }
   return maxint;
}

int zerodate = date2int("2020-02-24");

void mzcr2mzcr(const string& horizon, bool elementary=false)
{
    int lastdate = date2int(horizon);
    int numdates = lastdate - zerodate + 1;

    counter ocounter;
    unsigned inconsistento = 0;

    unsigned firstschool = 3;
    unsigned numclasses = 3+9+4+1;

    vector<vector<unsigned>> I(numdates,vector<unsigned>(numcohorts,0));
    vector<vector<unsigned>> R(numdates,vector<unsigned>(numcohorts,0));
    vector<vector<unsigned>> E(numdates,vector<unsigned>(numclasses,0));

    csv<','> osoby("/home/martin/Documents/s/covid/data/mzcr/osoby.csv");

    cout << "Importing osoby" << endl;

    for(unsigned i=1; i<osoby.r(); i++)
    {
        enum {datum,vek,pohlavi,kraj_nuts_kod,okres_lau_kod, nakaza_v_zahranici,nakaza_zeme_csu_kod};

        string ds = osoby(i,datum);
        int d = date2int(ds, zerodate,lastdate, ocounter);

        if(d == maxint)
            inconsistento++;
        else
        {
            const unsigned noage = 1000;
            unsigned v = noage;
            try
            {
                v = osoby.getunsigned(i,vek);
            }
            catch(...)
            {
                clog << "Bad age " << osoby(i,vek) << endl;
                inconsistento++;
            }
            if(v != noage)
            {
                cohorts c = v2cohort(v);

                if(v>=firstschool && v < firstschool + numclasses)
                    E[d-zerodate][v-firstschool]++;

                if(c==numcohorts)
                    inconsistento++;
                else
                {
                    R[d-zerodate][c]++;
                    if(osoby(i,nakaza_v_zahranici)=="1")
                        I[d-zerodate][c]++;
                }
            }
        }
    }
    cout << inconsistento << " records: "
         << ocounter.wrong << " wrong dates,"
         << ocounter.over << " dates over,"
         << ocounter.under << " dates under" << endl;

    counter ucounter;
    unsigned inconsistentu = 0;
    vector<vector<unsigned>> D(numdates,vector<unsigned>(numcohorts,0));

    csv<','> umrti("/home/martin/Documents/s/covid/data/mzcr/umrti.csv");

    cout << "Importing umrti" << endl;

    for(unsigned i=1; i<umrti.r(); i++)
    {
        enum labels {datum,vek,pohlavi,kraj_nuts_kod,okres_lau_kod};

        string ds = umrti(i,datum);
        int d = date2int(ds, zerodate,lastdate, ucounter);

        if(d == maxint)
            inconsistentu++;
        else
        {
            unsigned v = umrti.getunsigned(i,vek);
            cohorts c = v2cohort(v);
            if(c==numcohorts)
                inconsistentu++;
            else
                D[d-zerodate][c]++;
        }
    }
    cout << inconsistentu << " records: "
         << ucounter.wrong << " wrong dates,"
         << ucounter.over << " dates over,"
         << ucounter.under << " dates under" << endl;

    if(elementary)
    {
        cout << "S1,S2,S3,1,2,3,4,5,6,7,8,9,G1,G2,G3,G4,G5" << endl;
        for(unsigned i=0; i<numdates; i++)
        {
            for(unsigned j=0; j<numclasses; j++)
                cout << E[i][j] << ",";
            cout << endl;
        }
    }
    else
    {
        cout << "I0,I20,I65,I80,R0,R20,R65,R80,D0,D20,D65,D80" << endl;

        for(unsigned i=0; i<numdates; i++)
        {
            cout << I[i][0] << "," << I[i][1] << ","
                 << I[i][2] << "," << I[i][3] << ",";
            cout << R[i][0] << "," << R[i][1] << ","
                 << R[i][2] << "," << R[i][3] << ",";
            cout << D[i][0] << "," << D[i][1] << ","
                 << D[i][2] << "," << D[i][3] << endl;
        }
    }
}




// durs mean numbers to death
void olduzis2uzis(const string& horizon, bool uziscsv=false, bool rates=false, bool durs=false, bool hospdurs =false)
{
    enum cats {c0, c20, c65, numcats };

    enum labels {okres,orp,pohlavi,vek_kat,datum_prvniho_priznaku,
                 datum_odberu, datum_vysledku, datum_hlaseni,
                 datum_izolace, zahajeni_hospitalizace,
                 ukonceni_hospitalizace,datum_vyleceni,datum_umrti,
                 numlabels};

    enum hcatc {cy, co, numhcats };

    enum gender {man, woman, other, numgenders};

    int lastdate = date2int(horizon);
    int numdates = lastdate - zerodate + 1;

    vector<counter> cntr(numlabels); // i know the first will be unused
    unsigned inconsistent = 0;

    struct staterec { int min; int whofirst; cats cat; };

    vector<vector<staterec>> fromisgivendet(numdates);
    vector<vector<staterec>> fromisd(numdates);
    vector<vector<staterec>> fromh(numdates);
    vector<vector<double>> delays(numdates);

    vector<vector<double>> h(numhcats,vector<double>(numdates,0));
    vector<vector<double>> r(numhcats,vector<double>(numdates,0));
    vector<vector<double>> dhosp(numhcats,vector<double>(numdates,0));
    vector<vector<double>> dother(numhcats,vector<double>(numdates,0));

    const unsigned numdurs = 100;
    vector<vector<unsigned>> durhist(2*3,vector<unsigned>(numdurs,0));
    vector<vector<unsigned>> hospdurhist(2*3,vector<unsigned>(numdurs,0));

    enum efromis { efideath, efihosp, efidet};

    csv<';'> src("/home/martin/Documents/s/covid/data/epidemie/modely_05_datumy.csv");

    assert(src.c(0)==numlabels);

    for(unsigned i=1; i<src.r(); i++)
    {
        string cs = src(i,vek_kat);
        cats c;
        if(cs == "0-19")
            c = c0;
        else if(cs == "20-64")
            c = c20;
        else if(cs == "65+")
            c = c65;
        else
        {
            inconsistent++;
            continue;
        }

        gender g;
        auto pohlstr = src(i,pohlavi);
        if(pohlstr == "M")
            g = man;
        else if(pohlstr == "Z")
            g = woman;
        else
            g = other;

        vector<int> res(numlabels);
        for(unsigned j=datum_prvniho_priznaku; j<numlabels; j++)
            res[j] = date2int(src(i,j),zerodate,lastdate,cntr[j]);
        int dpp = res[datum_prvniho_priznaku];
        int dod = res[datum_odberu];
        int dh = res[zahajeni_hospitalizace];
        int dr = res[ukonceni_hospitalizace];
        int dd = res[datum_umrti];
        int drep = res[datum_hlaseni];

        int old = c == c65;
        if(dh < maxint)
            h[old][dh-zerodate]++;
        if(dr < maxint)
            r[old][dr-zerodate]++;
        if(dd < maxint)
        {
            if(dh < maxint)
                dhosp[old][dd-zerodate]++;
            else
                dother[old][dd-zerodate]++;
        }

        if(dod==maxint || dpp==maxint || dh < dpp)
            inconsistent++;
        else
        {
            if(dd < maxint && dd > dpp)
            {
                unsigned cindex = c;
                if(g == woman)
                    cindex += numcats;
                unsigned d = dd - dpp;
                if(d >= numdurs)
                    d = numdurs - 1;
                durhist[cindex][d]++;
            }

            if(dr < maxint && dr > dh)
            {
                unsigned cindex = c;
                if(g == woman)
                    cindex += numcats;
                unsigned d = dr - dh;
                if(d >= numdurs)
                    d = numdurs - 1;
                hospdurhist[cindex][d]++;
            }

            if(drep < maxint)
                delays[dpp-zerodate].push_back(drep-dpp);

            if(dpp<dod || (dpp==dod && dh<maxint)) // symp
            {
                int h = dh - dpp;
                int t = dod - dpp;
                int d = dd - dpp;
                int m = min(min(h,t),d);
                assert(m < maxint / 2);
                staterec r = { m,m==d ? 2 :(m==h ? 1 : 0),c};
                fromisgivendet[dpp-zerodate].push_back(r);

                if(dod < dh && dod < dd)
                {
                    int hdo = dh - dod;
                    int ddo = dd - dod;
                    m = min(hdo,ddo);
                    if(m < maxint /2)
                        fromisd[dod-zerodate].push_back({m, ddo < hdo ? 1 : 0,c});
                }
                if(dh < maxint)
                {
                    int d = dd - dh;
                    int r = dr - dh;
                    m = min(d,r);
                    if(m < maxint /2)
                        fromh[dh-zerodate].push_back({m, d < r ? 1 : 0,c});
                }
            }
        }
    }

    for(unsigned j=datum_prvniho_priznaku; j<numlabels; j++)
    {
        cout << src(0,j) << ": " << cntr[j].wrong << " wrong,"
               << cntr[j].over << " over,"
               << cntr[j].under << " under" << endl;
    }

    cout << inconsistent << " inconsistent records." << endl;
    cout << endl;

    if(durs)
    {
        for(unsigned i=0; i<numcats*2; i++)
        {
           cout << i;
           for(unsigned j=0; j<numdurs; j++)
               cout << "," << durhist[i][j];
           cout << endl;
        }
    }

    if(hospdurs)
    {
        for(unsigned i=0; i<numcats*2; i++)
        {
           cout << i;
           for(unsigned j=0; j<numdurs; j++)
               cout << "," << hospdurhist[i][j];
           cout << endl;
        }
    }


    if(uziscsv)
    {
        cout << "HY,HO,RHY,RHO,DHY,DHO,DOY,DOO" << endl;
        for(unsigned i=0; i<numdates; i++)
            cout
                 << h[0][i] << "," << h[1][i] << ","
                 << r[0][i] << "," << r[1][i] << ","
                 << dhosp[0][i] << "," << dhosp[1][i] << ","
                 << dother[0][i] << "," << dother[1][i] << endl;
    }

    if(rates)
    {
        cout << "no,";
        vector<string> cl = { "0", "20", "65", "all"};
        for(unsigned i=0; i<=numcats; i++)
        {
           cout << cl[i]+"_meanisgivendet, test, hosp, death,";
           cout << cl[i]+"_meanisd, hosp, death,";
           cout << cl[i]+"_meanh, recovery, death,";
        }
        cout << "repdelay" << endl;

        for(unsigned i=0; i < fromisgivendet.size(); i++ )
        {
            cout << i << ",";
            for(unsigned c=0; c <= numcats; c++)
            {
                bool all = c==numcats;
                double s=0;
                vector<double> ss(3,0);
                unsigned nall=fromisgivendet[i].size();
                unsigned n=0;
                for(unsigned j=0; j<nall; j++)
                {
                    const auto& r = fromisgivendet[i][j];
                    if(all || r.cat == c)
                    {
                        s += r.min;
                        ss[r.whofirst]++;
                        n++;
                    }
                }

                double sd=0;
                vector<double> ssd(2,0);
                unsigned ndall=fromisd[i].size();
                unsigned nd = 0;
                for(unsigned j=0; j<ndall; j++)
                {
                    const auto& r = fromisd[i][j];
                    if(all || r.cat == c)
                    {
                        sd += r.min;
                        ssd[r.whofirst]++;
                        nd++;
                    }
                }

                double sh=0;
                vector<double> ssh(2,0);
                unsigned nhall=fromh[i].size();
                unsigned nh = 0;
                for(unsigned j=0; j<nhall; j++)
                {
                    const auto& r = fromh[i][j];
                    if(all || r.cat == c)
                    {
                        sh += r.min;
                        ssh[r.whofirst]++;
                        nh++;
                    }
                }

                if(n)
                {
                    cout << s / n << ","<< ss[0] / n << ","
                         << ss[1] / n << "," << ss[2] / n << ",";
                }
                else
                    cout << ",,,,";

                if(nd)
                    cout << sd / nd << ","<< ssd[0] / nd << ","
                     << ssd[1] / nd << ",";
                else
                    cout << ",,,";

                if(nh)
                    cout << sh / nh << ","<< ssh[0] / nh << ","
                     << ssh[1] / nh << ",";
                else
                    cout << ",,,";

             }
             if(delays[i].size())
             {
                 double s=0;
                 for(unsigned j=0; j<delays[i].size(); j++)
                    s += delays[i][j];
                 cout << s /delays[i].size();
             }

             cout << endl;
        }
    }
}

void uzis2uzis(const string& horizon)
{


    enum cats {c0, c20, c65, c80, numcats };

    enum labels {vek_kat,pohlavi,kraj_bydliste,kraj_prvni_nemocnice,
                 datum_priznaku,datum_odberu,datum_positivity,stav_dle_khs,
                 zahajeni_hosp,
                 posledni_zaznam,stav_posledni_zaznam,posledni_hosp_zaznam,nejtezsi_stav,
                 tezky_stav,tezky_stav_pocatek,dni_tezky_stav,tezky_stav_posledni,
                 jip,jip_pocatek,dni_jip,jip_posledni,
                 kyslik,kyslik_pocatek,dni_kyslik,kyslik_posledni,
                 upv,upv_pocatek,dni_upv,upv_posledni,
                 ecmo,ecmo_pocatek,dni_ecmo,ecmo_posledni,
                 umrti,datum_umrti,numlabels};


    enum gender {man, woman, other, numgenders};

    enum dateinds {hosp, release, death, numdateinds};

    int lastdate = date2int(horizon);
    int numdates = lastdate - zerodate + 1;

    vector<counter> cntr(numdateinds); // i know the first will be unused
    unsigned inconsistent = 0;

/*    struct staterec { int min; int whofirst; cats cat; };

    vector<vector<staterec>> fromisgivendet(numdates);
    vector<vector<staterec>> fromisd(numdates);
    vector<vector<staterec>> fromh(numdates);
    vector<vector<double>> delays(numdates);*/

    vector<vector<double>> h(numcats,vector<double>(numdates,0));
    vector<vector<double>> r(numcats,vector<double>(numdates,0));
    vector<vector<double>> dhosp(numcats,vector<double>(numdates,0));

//    const unsigned numdurs = 100;
//    vector<vector<unsigned>> durhist(2*3,vector<unsigned>(numdurs,0));
//    vector<vector<unsigned>> hospdurhist(2*3,vector<unsigned>(numdurs,0));

//    enum efromis { efideath, efihosp, efidet};

    csv<';'> src("/home/martin/Documents/s/covid/data/epidemie/modely_05_hospitalizovani_analyza.csv");

    assert(src.c(0)==numlabels);

    for(unsigned i=1; i<src.r(); i++)
    {
        string cs = src(i,vek_kat);
        cats c;
        string ls;
        unsigned age;
        for(unsigned j=0; j<cs.size() && cs[j] != '-'; j++)
            ls += cs[j];
        try
        {
            age = stoul(ls);
        }
        catch (...)
        {
            clog << "Cannot convert'" << ls << "' to unsigned" << endl;
            inconsistent++;
            continue;

        }
        if(age >= 80)
            c = c80;
        else if(age >= 65)
            c = c65;
        else if(age >= 20)
            c = c20;
        else
            c = c0;

        gender g;
        auto pohlstr = src(i,pohlavi);
        if(pohlstr == "M")
            g = man;
        else if(pohlstr == "Z")
            g = woman;
        else
            g = other;

        vector<int> res(numlabels);
        int dh = date2int(src(i,zahajeni_hosp),zerodate,lastdate,cntr[hosp]);
        int dlr = date2int(src(i,posledni_zaznam),zerodate,lastdate,cntr[release]);
        int dd = date2int(src(i,datum_umrti),zerodate,lastdate,cntr[death]);

        if(dh < maxint)
            h[c][dh-zerodate]++;
        if(dd < maxint)
            dhosp[c][dd-zerodate]++;
        else if(dlr < maxint && src(i,posledni_hosp_zaznam)[0] != 'h') // hospitalizace okracuje
            r[c][dlr-zerodate]++;

        if(dh==maxint)
        {
            clog << "record without hosp date" << endl;
            inconsistent++;
        }
        else if(dlr==maxint)
        {
            clog << "record without last record or " << endl;
            inconsistent++;
        }
    }

    for(unsigned j=0; j<numdateinds; j++)
    {
        cout << j << ": " << cntr[j].wrong << " wrong,"
               << cntr[j].over << " over,"
               << cntr[j].under << " under" << endl;
    }

    cout << inconsistent << " inconsistent records." << endl;
    cout << endl;


    cout << "H0,H20,H65,H80,R0,R20,R65,R80,DH0,DH20,DH65,DH80" << endl;
    for(int i=0; i<numdates; i++)
    {
        for(unsigned j=0; j<numcats; j++)
            cout << h[j][i] << ",";
        for(unsigned j=0; j<numcats; j++)
            cout << r[j][i] << ",";
        for(unsigned j=0; j<numcats; j++)
            cout << dhosp[j][i] << ",";
        cout << endl;
    }

}


void mzcr2districts(const string& horizon)
{
    int lastdate = date2int(horizon);
    int numdates = lastdate - zerodate + 1;
//    int firstwnum = 9;
    int numweeks = (numdates + 6) / 7;
    int lastweekinoutput = 72;
    if(lastweekinoutput < numweeks)
        throw "too much weeks";
    counter ocounter;
    unsigned inconsistento = 0;

    unsigned firstschool = 3;
    unsigned numclasses = 3+9+4+1;
    unsigned numdistricts = districts.size();


    vector<vector<unsigned>> C(numdistricts,vector<unsigned>(numweeks,0));
    vector<vector<vector<unsigned>>> E(numdistricts,
             vector<vector<unsigned>>(numdates,vector<unsigned>(numclasses,0)));
    vector<unsigned> Y(numweeks,0);
    vector<unsigned> YM(numweeks,0);

    csv<','> osoby("/home/martin/Documents/s/covid/data/mzcr/osoby.csv");

    cout << "Importing osoby" << endl;

    for(unsigned i=1; i<osoby.r(); i++)
    {
        enum {datum,vek,pohlavi,kraj_nuts_kod,okres_lau_kod, nakaza_v_zahranici,nakaza_zeme_csu_kod};

        string ds = osoby(i,datum);
        int d = date2int(ds, zerodate,lastdate, ocounter);

        if(d == maxint)
        {
            clog << "Invalid or out of range date " << ds << endl;
            inconsistento++;
        }
        else
        {
            int w = (d-zerodate) / 7;
            const unsigned noage = 1000;
            unsigned v = noage;
            try
            {
                v = osoby.getunsigned(i,vek);
            }
            catch(...)
            {
                clog << "Bad age " << osoby(i,vek) << endl;
                inconsistento++;
            }
            if(v != noage)
            {
                int distr = str2district(osoby(i,okres_lau_kod));
                if(distr == maxint)
                {
                    clog << "Unknown district code " << osoby(i,okres_lau_kod) << endl;
                    inconsistento++;
                }
                else
                {
                    if(v>=firstschool && v < firstschool + numclasses)
                        E[distr][w][v-firstschool]++;
                    else
                        YM[w]++;
                    C[distr][w]++;
                    Y[w]++;
                }
            }
        }
    }
    cout << inconsistento << " records: "
         << ocounter.wrong << " wrong dates,"
         << ocounter.over << " dates over,"
         << ocounter.under << " dates under" << endl;

    ofstream o(sys::outputfolder()+"schools.csv");
    if(!o)
    {
        cerr << "Cannot open " + sys::outputfolder()+"schools.csv" << endl;
        throw;
    }

    o << "D,W,Y,YM,YD,I,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19" << endl;
    for(unsigned i=0; i<numdistricts; i++)
    {
        int j=1;
        for(; j<numweeks; j++) // we start on 2.3.
        {
            o << i << "," << j << "," << Y[j] << "," << YM[j] << "," << C[i][j] << "," << districts[i].inhabitants;
            for(unsigned k=0; k<numclasses; k++)
                o << "," << E[i][j][k];
            o << endl;
        }
        for(; j<=lastweekinoutput; j++)
            o << i << "," << j << endl;
    }
}

