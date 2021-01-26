#include "orpp/csv.hpp"
#include <limits.h>
using namespace std;
using namespace orpp;

constexpr int maxint = numeric_limits<int>::max();


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


void uzis2uzis()
{
    enum labels {okres,orp,pohlavi,vek_kat,datum_prvniho_priznaku,
                 datum_odberu, datum_vysledku, datum_hlaseni,
                 datum_izolace, zahajeni_hospitalizace,
                 ukonceni_hospitalizace,datum_vyleceni,datum_umrti,
                 numlabels};



    enum cats {c0, c20, c65, numcats };
    enum hcatc {cy, co, numhcats };

    enum gender {man, woman, other, numgenders};


    int zerodate = date2int("2020-02-24");
    int lastdate = date2int("2021-01-16");
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

    enum efromis { efideath, efihosp, efidet};

    csv<';'> src("/home/martin/tmp/pepa/IDEA-anti-COVID-19-data/epidemie/modely_05_datumy.csv");

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

    if(0)
    {
        for(unsigned i=0; i<numcats*2; i++)
        {
           cout << i;
           for(unsigned j=0; j<numdurs; j++)
               cout << "," << durhist[i][j];
           cout << endl;
        }
    }

    if(1)
    {
        cout << "HY,HO,RHY,RHO,DHY,DHO,DOY,DOO" << endl;
        for(unsigned i=0; i<numdates; i++)
            cout
                 << h[0][i] << "," << h[1][i] << ","
                 << r[0][i] << "," << r[1][i] << ","
                 << dhosp[0][i] << "," << dhosp[1][i] << ","
                 << dother[0][i] << "," << dother[1][i] << endl;
    }

    if(0)
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
