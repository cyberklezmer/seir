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
    vector<vector<double>> data;
    unsigned Deltad=10;
    unsigned trsh=10;
    double c0 = 2000;

    vector<double> z;
    vector<vector<double>> reg;
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
            for(unsigned k=4; k<input[0].c(i); k++)
            {
                vector<double> rec;
                bool isnz = false;
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
                        clog << "first nonzero " << k-4 << endl;
                        isnz = true;
                    }
                }
                if(isnz)
                    c.data.push_back(rec);
            }
            for(unsigned i=0;i<c.data.size(); i++)
            {
                 c.data[i][eI] -= c.data[i][eR] + c.data[i][eD];
                 for(unsigned j=0; j<c.data[i].size(); j++)
                      clog << c.data[i][j] << ",";
                 clog << endl;
            }
            return;
        }
    }
    es  << " not found!" << endl;
    throw(es.str());
}

enum eparams { ealpha, edeltaa, egammaa, egammai,ephii, epsii, emu1, etheta1, emu2, etheta2, enumpars };

const vector<string> parnames = {"alpha", "deltaa","gammaa","gammai", "phii", "psii",  "mu1", "theta1","mu2","theta2"};
const vector<double> initvals = {20, 0.2,0.1,0.1, 0.01, 0.05, 0.5, 0.5, 0.5, 0.5};


double olsobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
{    
    country& c = *(static_cast<country*>(f_data));
    c.reg.resize(0);
    c.z.resize(0);
    const vector<vector<double>> x = c.data;
    double r = 0.0;
    vector<double> last(enumcols,0.0);
    double ztm1 = v[ealpha];
//for(unsigned i=0;i<enumpars; i++)
//    clog << v[i] << ", v=";
//clog << endl;
    for(unsigned t=0; t<x.size(); t++)
    {
        double ytm1;
        if(t >=c.Deltad+1)
        {
            double Itau = x[t-(c.Deltad+1)][eI];
            ytm1 = Itau * v[ephii] + max(Itau - c.c0,0.0) * v[epsii];
        }
        else
            ytm1 = 0;
        double thetatm1 = t < c.trsh+1 ? v[etheta1] : v[etheta2];
        double mutm1 = t < c.trsh+1 ? v[emu1] : v[emu2];
        double nutm1 = 1-v[edeltaa]-v[egammaa] + thetatm1;
        if(t>=2)
            ztm1 += nutm1*ztm1 + mutm1 * last[eI];
        const vector<double> pres = x[t];
        vector<double> reg =
        {
            (1-v[egammai]) * last[eI] + v[edeltaa] * ztm1 - ytm1,
            last[eD]+ ytm1 ,
            last[eR] + v[egammaa] * ztm1 + v[egammai] * last[eI]
        };

        double dI = reg[eI] - pres[eI];
        double dD = reg[eD] - pres[eD];
        double dR = reg[eR] - pres[eR];

        double wI = max(fabs(v[egammai] * last[eI] + v[edeltaa] * ztm1),10.0);
        double wD = max(ytm1,10.0);
        double wR = max(fabs(v[egammaa] * ztm1 + v[egammai] * last[eI]),10.0);

        r += (dI*dI) / wI + dD*dD / wD +dR*dR / wR;
        last = pres;
        c.reg.push_back(reg);
        c.z.push_back(ztm1);
    }
//clog << r << endl;
    return r;
}

vector<double> ols(const country& ac, vector<paraminfo> p)
{
    country c(ac);
    opt o(LN_COBYLA,enumpars);
    o.set_lower_bounds(vector<double>(enumpars,0));
    o.set_ftol_abs(0.00001);
    o.set_maxtime(10);

    o.set_min_objective(olsobj, &c);

    vector<double> v(enumpars);
    for(unsigned i=0; i<enumpars; i++)
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
    for(unsigned i=0; i<enumpars; i++)
    {
        cout << p[i].name << "=" << v[i] << endl;
    }

    assert(c.data.size()==c.reg.size());
    assert(c.data.size()==c.z.size());
    for(unsigned i=0; i<c.data.size(); i++)
    {
        clog << c.data[i][eI] << "," << c.reg[i][eI] << ","
             << c.data[i][eD] << "," << c.reg[i][eD] << ","
             << c.data[i][eR] << "," << c.reg[i][eR] << ","
             << c.z[i] << endl;
    }
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

        vector<paraminfo> p(enumpars);
        for(unsigned i=0; i<enumpars; i++)
        {
            p[i].name = parnames[i];
            p[i].initial = initvals[i];
        }
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
