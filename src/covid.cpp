#include "cctools.hpp"
#include "epp.hpp"

using namespace std;
using namespace cctools;
using namespace epp;
using namespace nlopt;


enum ecolumn { eI, eD, eR, enumcols};

const string datafolder = "../data/COVID-19/csse_covid_19_data/csse_covid_19_time_series/";

const vector<string> csvn =
{
  "time_series_covid19_confirmed_global.csv",
  "time_series_covid19_deaths_global.csv",
  "time_series_covid19_recovered_global.csv"
};

const string countriesfn = "../src/countries.csv";

struct country
{
    static const unsigned jan22 = 43852;

    country(const string& an, const string& ap) :
        province(ap), name(an) {}
    string province;
    string name;
    vector<vector<double>> data;
    vector<double> a;
    vector<vector<double>> reg;
    unsigned datastart;
    unsigned lockdown;
    unsigned nums() const
    {
       return lockdown < datastart ? 0 : lockdown-datastart;
    }
};


country addcountry(const vector<csv<','>>& input, const string& name,
                const string& province)
{
    country c(name,province);
    ostringstream es;
    es << "Country " << c.name << " (" << c.province << "): ";

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
        csv<','> cf(countriesfn);
        if(cf.r()<2 || cf.c(0)<3 )
            throw "stg wrong with countries";
clog << cf(i,1) << endl;
        if(cf(i,0)==province && cf(i,1)==name)
        {
            unsigned long v;
            if(cf.c(i)>=3)
            {
               try
               {
                   v = stoul(cf(i,2));
               }
               catch (...)
               {
                   es << " col " << i << " cannot convert " << cf(i,2);
                   throw es.str();
               }
               c.lockdown = v;
            }
            else
                c.lockdown = country::jan22;
        }
    }
    assert(c.data.size()==0);
    vector<vector<double>> data;
    clog << "Adding " << c.name << " (" << c.province << ")" << endl;

    vector<unsigned> indices;
    for(unsigned k=0;k<3;k++)
    {
        bool found = false;
        for(unsigned i=1; i<input[1].r(); i++)
        {
            if(input[k](i,0) == c.province && input[k](i,1) == c.name)
            {
                indices.push_back(i);
                found = true;
                break;
            }
        }
        if(!found)
        {
            es << "cannot find for k=" << k;
            throw(es.str());
        }
    }

    assert(indices.size()==3);
    unsigned t=1;
    c.datastart = country::jan22;
    bool isnz = false;
    for(unsigned k=4; k<input[0].c(indices[0]); k++, t++)
    {
        vector<double> rec;
        for(unsigned j=0; j<input.size(); j++)
        {
            unsigned i = indices[j];
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
                clog << "first nonzero " << c.datastart << endl;
                isnz = true;
            }
        }
        if(isnz)
{
clog << rec[0] << endl;
            data.push_back(rec);
}
        else
            c.datastart++;
    }

    c.data.push_back(data[0]);
    for(unsigned i=1;i<data.size(); i++)
    {
         vector<double> rec;
         for(unsigned j=0; j<data[i].size(); j++)
         {
              rec.push_back(data[i][j]-data[i-1][j]);
              clog << rec[j] << ",";
         }
         clog << endl;
         c.data.push_back(rec);
    }
    clog << "Added with s=" << c.nums() << endl;
    return c;
}

enum eparams { eiotaalpha, etheta, etheta2, eiotanu, enumallpars };


const vector<string> parnames =
{"iotaalpha", "theta", "theta2", "iotanu" };


const vector<double> upper =
{ 100,        3,       3,       100 };

class covidmle : public mle
{
    const country& c;
vector<double> yy;

    vector<double> Ik;
    vector<double> y;
    vector<double> thetader;
    vector<double> theta1der;
    vector<double> iotaalphader;
    vector<double> iotanuder;


    virtual int getN()
    {
        return c.data.size();
    }
public:

    virtual void beforeloglikeval(const vector<double> ap)
    {
//clog << "(" << ap[0] << "," << ap[1] << ")" << endl;
/*
yy.resize(0);
vector<double> p;
for(unsigned i=0; i<enumallpars; i++)
    if(i < ap.size())
        p.push_back(ap[i]);
    else
        p.push_back(0);
        // all starting with paper's zero
        const vector<vector<double>> x = c.data;

        iotaalphader.resize(0);
        iotanuder.resize(0);
        thetader.resize(0);
        theta1der.resize(0);
        Ik.resize(0);
        y.resize(0);

        y.push_back(p[eiotaalpha]);
yy.push_back(p[eiotaalpha]);
        Ik.push_back(0);
        thetader.push_back(0);
        theta1der.push_back(0);
        iotaalphader.push_back(1);
        iotanuder.push_back(0);

        for(unsigned t=1; t<x.size(); t++)
        {
            int tau = t-1;
            double I=0;
            for(unsigned i=0; i<c.k; i++,tau--)
            {
                if(tau >= 0)
                    I += x[tau][eI];
            }
            Ik.push_back(I);

            bool withs = t <= c.s;

            if(withs)
            {
                y.push_back(p[etheta]*y[t-1] + p[eiotanu] * Ik[t-1]);
                thetader.push_back(y[t-1] + p[etheta] * thetader[t-1]);
                theta1der.push_back(0);
                iotaalphader.push_back(iotaalphader[t-1] *p[etheta]);
yy.push_back((p[etheta])*yy[t-1] + p[eiotanu] * Ik[t-1]);
                iotanuder.push_back(p[etheta] *iotanuder[t-1] + Ik[t-1]);
            }
            else
            {
                y.push_back(p[etheta2]*y[t-1] + p[eiotanu] * Ik[t-1]);
yy.push_back((p[etheta2]+0.0001)*yy[t-1] + p[eiotanu] * Ik[t-1]);
                thetader.push_back( p[etheta2] * thetader[t-1]);
                theta1der.push_back(y[t-1] + p[etheta2] * theta1der[t-1]);
                iotaalphader.push_back(iotaalphader[t-1] *p[etheta2]);
                iotanuder.push_back(p[etheta2] *iotanuder[t-1] + Ik[t-1]);
*/
//if(t==60)
//  clog << "!!!!"  << (yy[yy.size()-1]-y[y.size()-1])/0.0001 << "=" << theta1der[iotaalphader.size()-1] << endl;
//            }
//clog << t << " y=" << y[y.size()-1] << " td= " << thetader[thetader.size()-1] <<  "..."
//      << Ik[t-1] <<  endl;
//        }
//throw;

    }

//    virtual void afterloglikeval(const vector<double>& aparams,
//                         const vector<double>& grad, double loglik)
//    {
//clog << "ll=" << loglik << ": ";
//for(unsigned i=0; i<grad.size(); i++)
//    clog << aparams[i] << "(" << grad[i] << "),";
//clog << endl;
//    }

    virtual double evallogdensity(int ii, const vector<double>& p,
                                                  vector<double>& g)
    {
        unsigned i = static_cast<unsigned int>(ii);  // to mute warnings
        const vector<vector<double>> x = c.data;
        unsigned s = c.nums();
        unsigned nums = min(i,s);
        unsigned numt = max(i,s)-s;
        double d = pow(p[etheta],nums) * pow(p[etheta2],numt);
        double z = p[eiotaalpha] * d;
        double iader = d;
        double thetader = nums ? p[eiotaalpha] * nums * pow(p[etheta],nums-1) * pow(p[etheta2],numt) : 0;
        double theta2der = numt ? p[eiotaalpha] * pow(p[etheta],nums) * numt * pow(p[etheta2],numt-1) : 0;
        double dif = z - x[i][eI];
        double w = max(fabs(x[i][eI]),1.0);

        g[eiotaalpha] = - 1.0 / w * dif * iader;
        g[etheta] = - 1.0 / w * dif * thetader;
        g[etheta2] = - 1.0 / w * dif * theta2der;

double dd =  pow(p[etheta],nums) * pow((p[etheta2]+0.0001),numt);
//clog <<  p[eiotaalpha]*(dd - d) /0.0001 << "=*=" << theta2der << endl;
double difdif = p[eiotaalpha] * dd - x[i][eI];
//        g[eiotanu] = - 1.0 / w * dif * iotanuder[i];
/*        g[etheta2] = - 1.0 / w * dif * theta1der[i];


        double dif = y[i] - x[i][eI];
        double w = max(fabs(y[i]),1.0);
double difdif = yy[i]-x[i][eI];
        g[eiotaalpha] = - 1.0 / w * dif * iotaalphader[i];
        g[etheta] = - 1.0 / w * dif * thetader[i];

//        g[eiotanu] = - 1.0 / w * dif * iotanuder[i];
        g[etheta2] = - 1.0 / w * dif * theta1der[i]; */
//clog << (- 1.0 / (2 * w) * difdif * difdif - (- 1.0 / (2 * w) * dif * dif) ) / 0.0001
//    << "=dd=" << g[etheta2] << endl;
        return - 1.0 / (2 * w) * dif * dif;
    }
public:


    covidmle(const vector<paraminfo>& aparams, const country& ac):
        c(ac), mle(aparams)
      {  }
};




const unsigned pred = 20;

double olsobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
{
    country& c = *(static_cast<country*>(f_data));
    c.reg.resize(0);
    c.a.resize(0);
    const vector<vector<double>> x = c.data;
    double r = 0.0;

    vector<double> last(enumcols,0.0);
    double iotaztm1 = v[eiotaalpha];

    for(unsigned t=0; t<x.size()+pred; t++)
    {
        vector<double> reg =
        {
            iotaztm1
        };

        if(t<x.size())
        {
            const vector<double> pres = x[t];

            double dI = reg[eI] - pres[eI];

            double wI = max(fabs(iotaztm1),5.0);

            r += (dI*dI) / wI;
            c.reg.push_back(reg);
        }
        else
        {
            c.reg.push_back(reg);
        }

        double lastIActive = 0;
        int tau=t-1;
        for(unsigned i=0; i< 0/*c.k*/; i++,tau--)
        {
            if(tau >= 0 && tau>=x.size())
                lastIActive += c.reg[tau][eI];
            else if(tau >= 0)
                lastIActive += x[tau][eI];
        }
        double thetatm1 = t+1 <= c.nums() ? v[etheta] : v[etheta2];
        iotaztm1 = thetatm1*iotaztm1 + v[eiotanu] * lastIActive ;
        if(t<x.size())
            last =  x[t];
        else
            last = reg;

        c.a.push_back(iotaztm1);
    }
    return r / 2;
}



vector<double> ols(const country& ac, vector<paraminfo> p)
{
    country c(ac);
    opt o(LN_COBYLA,enumallpars);
    o.set_lower_bounds(vector<double>(enumallpars,0));
    o.set_upper_bounds(upper);
    o.set_ftol_abs(1e-10);
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
    assert(c.data.size()+pred==c.a.size());
    unsigned i=0;
    for(; i<c.data.size(); i++)
    {
        clog << c.data[i][eI] << "," << c.reg[i][eI] << ","
             << c.data[i][eD] << "," << c.reg[i][eD] << ","
             << c.data[i][eR] << "," << c.reg[i][eR] << ","
             << c.a[i] << endl;
    }
    for(;i<c.data.size()+pred; i++)
    {
        clog << "," << c.reg[i][eI] << ","
             << "," << c.reg[i][eD] << ","
             << "," << c.reg[i][eR] << ","
             << c.a[i] << endl;
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
//        country c =addcountry(input, "Italy","");
//        country c =addcountry(input, "China","Hubei");
        country c =addcountry(input, "Korea,South","");

        vector<paraminfo> p(enumallpars);
        for(unsigned i=0; i<enumallpars; i++)
        {
            p[i].name = parnames[i];
        }

        double mu = 0.2;
        double nu = 0.1;
        double iota = 0.2;
        double gamma = 1.0 / 12.0;
        double delta = 0.03;
        double gammai = 0.11;
        double phi = 1.0 / 15.0;
        double gammac = 2.0 / 15.0;
        double alpha = 1;

//        p[esigma].initial = 1;
        p[eiotaalpha].initial = iota * alpha ;
        p[etheta].initial= 1+mu-gamma-iota;
        p[etheta2].initial= 1+mu-gamma-iota;
        p[eiotanu].initial= iota*nu;

//p[eiotaalpha].initial = 1 ;
//p[etheta].initial= 0.5;
//p[etheta2].initial= 0.8;
//p[eiotanu].initial=0.2;


//        p[ephidelta].initial= phi*delta;
//        p[eomega].initial=1 - delta - gammac;
//        p[egammai].initial= gammai;
//        p[egammacdelta].initial= gammac * delta;

        vector<double> olsres = ols(c,p);

        vector<paramresult> dr;
        double ov = -HUGE_VAL;

        p.resize(3);
        double net = 3;
        vector<vector<double>> J;
        for(unsigned i=1; i<net; i++)
            for(unsigned j=1; j<net; j++)
            {
                for(unsigned i=0; i<p.size(); i++)
                {

                    p[i].lower = 0;
                    p[i].upper = upper[i];
                    p[i].initial = olsres[i];
                }

                p[etheta].initial = i*p[etheta].upper / net;
                p[etheta2].initial = j*p[etheta2].upper / net;
                vector<paramresult> r;
                covidmle m(p,c);
                m.setftolabs(1e-6);
                double res = m.estimate(r,true);
                clog << res << endl;
                if(res > ov)
                {
                    ov = res;
                    dr = r;
                    J = m.getJ();
                }
            }

        clog << "ov=" << ov << endl;
        for(unsigned i=0; i<dr.size(); i++)
        {
            clog << dr[i] << endl;
        }

        ostringstream fn;
        fn << c.name << c.province << ".csv";
        ofstream out(fn.str());
        if(!out)
        {
            es << "Cannot open '" << fn.str() << "'";
            throw es.str();
        }
        out << "T,dAe,dI,dIe,dD,dDe,dR,dRe" << endl;
        double iotaa = dr[eiotaalpha].value;
        vector<double> last(3,0.0);

        for(unsigned i=0; i<c.data.size()+pred; i++)
        {
            bool pred = i>=c.data.size();
            unsigned t = c.datastart+i;
            double lastiotaa = iotaa;
            iotaa *= t <= c.lockdown ? dr[etheta].value : dr[etheta2].value;
            out << t << "," << iotaa << ",";
            if(!pred)
                out << c.data[i][eI] << "," << lastiotaa << ","
                    << c.data[i][eD] << ",,"
                    << c.data[i][eR] << ",," ;
            else
                out << ",," << lastiotaa << ",,,,," ;
            if(t==c.lockdown)
                out << "L";
            out << endl;
        }



/*        vector<double> xx = { dr[0].value, dr[1].value, dr[2].value};


        throw;
covidmle m(p,countries[0]);
double ll, oldll;

for(double d = -0.001 ; d < 0.001; d+= 0.0001)
{
    vector<double> g(3);
    xx[0]=dr[0].value + d;

    ll=m.loglik(xx,g);
    clog << "ll(" << xx[0]<<","<< xx[1]<<"," << xx[2] << ") "<< "=" << ll << ", g=," << g[2] << ", " << endl;
    if(d>-0.001+0.00000001)
        clog << oldll << "-" << ll
               << "..." << (oldll-ll)/0.0001 << ",=covid=" << g[0] << endl;
    else
        clog << endl;
    oldll = ll;
}*/
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


double sample[] = {1,2,5,3,2,6,1,8,1};

class testmle : public mle
{
    virtual int getN()
    {
        return sizeof(sample)/sizeof(sample[0]);
    }
public:
    virtual void afterloglikeval(const vector<double>& aparams,
                         const vector<double>& grad, double loglik)
    {
clog << grad[0] << "," << grad[1] << endl;

    }


    virtual double evallogdensity(int i, const vector<double>& x,
                                                  vector<double>& g)
    {
        double s = (sample[i]-x[0])*(sample[i]-x[0]);
        g[0] = (sample[i]-x[0]) / x[1] ;
        g[1] = - 1.0 / x[1] / 2 + s / (x[1] * x[1])/2.0;
        return -log(2*3.1415*x[1]) / 2.0 - s / x[1] / 2.0;
    }
public:

    testmle(const vector<paraminfo>& aparams):
        mle(aparams)
      {}
};


int _main()
{
    vector<paraminfo> p(2);
    p[0].name = "m";
    p[0].initial = 1;
    p[1].name = "V";
    p[1].lower = 0.00001;
    p[1].initial = 1;

    testmle m(p);
    m.setxtolabs(1e-10);

    vector<paramresult> r;

    m.estimate(r,true);
    cout.precision(15);
    cout << r[0] << endl << r[1] << endl;
}
