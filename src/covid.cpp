#include "cctools.hpp"
#include "epp.hpp"

using namespace std;
using namespace cctools;
using namespace epp;
using namespace nlopt;

enum ecolumn { eI, eD, eR, enumcols};

const string datafolder = "../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/";

const vector<string> csvn =
{ "time_series_19-covid-Confirmed.csv",
  "time_series_19-covid-Deaths.csv",
  "time_series_19-covid-Recovered.csv"
};

const string countriesfn = "../src/countries.csv";

struct country
{
    country(const string& ap, const string& an) :
        province(ap), name(an) {}
    string province;
    string name;
    unsigned firstnz = 1;
    vector<vector<double>> data;
    vector<double> cr;
    vector<double> cd;
    vector<double> a;
    vector<vector<double>> reg;
    unsigned s = 10;
};




void addcountry(const vector<csv<','>>& input, country& c)
{
    ostringstream es;
    es << "Country " << c.name << " (" << c.province << "): ";
    assert(c.data.size()==0);

    for(unsigned i=1; i<input[1].r(); i++)
    {
        if(input[0](i,0) == c.province && input[0](i,1) == c.name)
        {
            for(unsigned j=0; j<input.size(); j++)
            {
                if(!(input[j](i,0) == c.province && input[0](i,1) == c.name))
                {
                    es << "inconsistent input";
                    throw(es.str());
                }
            }
            unsigned t=1;
            bool isnz = false;
            for(unsigned k=4; k<input[0].c(i); k++, t++)
            {
                vector<double> rec;
                for(unsigned j=0; j<input.size(); j++)
                {
                    if(input[j].c(i)<k)
                        throw "different length of datasets";
                    unsigned long v;
                    try
                    {
                        v = stoul(input[j](i,k));
                    }
                    catch (...)
                    {
                        es << " col " << k << " cannot convert " << input[j](i,k);
                        throw es.str();
                    }                    
                    rec.push_back(v);
                    if(v>0 && !isnz)
                    {
                        c.firstnz = t;
                        clog << "first nonzero " << t << endl;
                        isnz = true;
                    }
                }
                if(isnz)
                    c.data.push_back(rec);
            }
            for(unsigned i=0;i<c.data.size(); i++)
            {
                 for(unsigned j=0; j<c.data[i].size(); j++)
                      clog << c.data[i][j] << ",";
                 clog << c.data[i][eI] - c.data[i][eR] - c.data[i][eD] << endl;
            }
            return;
        }
    }
    es  << " not found!" << endl;
    throw(es.str());
}

enum eparams { esigma, eiotaalpha, etheta, etheta2, eiotanu, ephidelta, eiparams=ephidelta, eomega, egammai , egammacdelta, enumallpars };


const vector<string> parnames =
{"sigma", "iotaalpha", "theta", "theta2", "iotanu",  "phidelta", "omega", "gammai" , "gammacdelta"};


const vector<double> upper =
{ 1,       1000,        2,       2,        1,         1,         1 ,       1,         1 };



const unsigned pred = 20;

double olsobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
{    
    country& c = *(static_cast<country*>(f_data));
    c.reg.resize(0);
    c.a.resize(0);
    c.cr.resize(0);
    c.cd.resize(0);
    const vector<vector<double>> x = c.data;
    double r = 0.0;

    vector<double> last(enumcols,0.0);
    double iotaztm1 = v[eiotaalpha];
    double phiytm1 = 0;
    double gammacytm1 = 0;

    for(unsigned t=0; t<x.size()+pred; t++)
    {
        vector<double> reg =
        {
            v[esigma] * last[eI] + iotaztm1,
            last[eD] + phiytm1,
            last[eR] + v[egammai] * last[eI] + gammacytm1
        };

        if(t<x.size())
        {
            const vector<double> pres = x[t];

            double dI = reg[eI] - pres[eI];
            double dD = reg[eD] - pres[eD];
            double dR = reg[eR] - pres[eR];

            double wI = max(fabs((v[esigma]-1)* last[eI]) + fabs(iotaztm1),5.0);
            double wD = max(fabs(phiytm1),1.0);
            double wR = max(fabs(v[egammai] * last[eI]) + fabs(gammacytm1),1.0);

            r += (dI*dI) / wI + (dD*dD) / wD + (dR*dR) / wR;
            c.reg.push_back(reg);
        }
        else
        {
            c.reg.push_back(reg);
        }

        double lastIActive = last[eI]-last[eR]-last[eD];
        double thetatm1 = t+1 <= c.s ? v[etheta] : v[etheta2];
        iotaztm1 = thetatm1*iotaztm1 + v[eiotanu] * lastIActive ;
        phiytm1 = v[eomega]*phiytm1 + v[ephidelta] * lastIActive;
        gammacytm1 = v[eomega]*gammacytm1 + v[egammacdelta] * lastIActive;

        if(t<x.size())
            last =  x[t];
        else
            last = reg;

        c.a.push_back(iotaztm1);
        c.cr.push_back(gammacytm1);
        c.cd.push_back(phiytm1);
    }
//clog << r << endl;
    return r;
}

vector<double> ols(const country& ac, vector<paraminfo> p)
{
    country c(ac);
    opt o(LN_COBYLA,enumallpars);
    o.set_lower_bounds(vector<double>(enumallpars,0));
    o.set_upper_bounds(upper);
    o.set_ftol_abs(0.00001);
    o.set_maxtime(10);

    o.set_min_objective(olsobj, &c);

    vector<double> v(enumallpars);
    for(unsigned i=0; i<enumallpars; i++)
        v[i] = p[i].initial;
    nlopt::result oresult;
    double ov;
    try
    {
        oresult = o.optimize(v, ov);
    }
    catch(std::exception &e)
    {
        throw e.what();
    }
    if(oresult < 0)
    {
        ostringstream es;
        es << "Negative result of nlopt: " << oresult << endl;
        throw es.str();
    }
    assert(c.data.size()+pred==c.reg.size());
    assert(c.data.size()+pred==c.cr.size());
    assert(c.data.size()+pred==c.cd.size());
    assert(c.data.size()+pred==c.a.size());
    unsigned i=0;
    for(; i<c.data.size(); i++)
    {
        clog << c.data[i][eI] << "," << c.reg[i][eI] << ","
             << c.data[i][eD] << "," << c.reg[i][eD] << ","
             << c.data[i][eR] << "," << c.reg[i][eR] << ","
             << c.a[i] << "," << c.cd[i] << "," << c.cr[i] << endl;
    }
    for(;i<c.data.size()+pred; i++)
    {
        clog << "," << c.reg[i][eI] << ","
             << "," << c.reg[i][eD] << ","
             << "," << c.reg[i][eR] << ","
             << c.a[i] << "," << c.cd[i] << "," << c.cr[i] << endl;
    }
    for(unsigned i=0; i<enumallpars; i++)
    {
        clog << p[i].name << "=" << v[i] << endl;
    }
    clog << "ov=" << ov << endl;
    return v;
}

int main(/* int argc, char *argv[] */)
{
    try
    {
        vector<csv<','>> input;
        for(unsigned i=0; i<enumcols; i++)
        {
            ostringstream s;
            s << datafolder << csvn[i];
            input.push_back(s.str());
        }

        ostringstream es;
        vector<country> countries;
        csv<','> cf(countriesfn);
        if(cf.r()<2 || cf.c(0)<3 )
            throw "stg wrong with countries";
        for(unsigned i=1; i<cf.r(); i++)
        {
            if(cf.c(i) < 2)
            {
                es << "Little cols at line " << i << endl;
                throw(es.str());
            }
            if(cf(i,2)=="1")
                countries.push_back(country(cf(i,0),cf(i,1)));
        }
        for(unsigned i=0; i<countries.size(); i++)
            addcountry(input, countries[i]);

        vector<paraminfo> p(enumallpars);
        for(unsigned i=0; i<enumallpars; i++)
        {
            p[i].name = parnames[i];
        }

        double mu = 1;
        double nu = 0.2;
        double iota = 0.2;
        double gamma = 1.0 / 12.0;
        double delta = 0.03;
        double gammai = 0.11;
        double phi = 1.0 / 15.0;
        double gammac = 2.0 / 15.0;
        double alpha = 100;

        p[esigma].initial = 1;
        p[eiotaalpha].initial = iota * alpha ;
        p[etheta].initial= 1+mu-gamma-iota;
        p[etheta2].initial= 1+mu-gamma-iota;
        p[eiotanu].initial= iota*nu;
        p[ephidelta].initial= phi*delta;
        p[eomega].initial=1 - delta - gammac;
        p[egammai].initial= gammai;
        p[egammacdelta].initial= gammac * delta;

        ols(countries[0],p);


    }
    catch (exception& e)
    {
         cerr << e.what() << endl;
     }
     catch (const char* m)
     {
            cerr << m << endl;
     }
     catch (const string& m)
     {
            cerr << m << endl;
     }

}
