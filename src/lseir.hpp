#ifndef LSEIR_HPP
#define LSEIR_HPP

#include "seirfilter.hpp"
#include "estimableseirfilter.hpp"
#include "orpp/csv.hpp"
#include <cmath>

class lseir: virtual public seirfilter
{
public:
    static constexpr unsigned fineobs = 40;

    enum states { E,  Ia,  Ip,  Is, numactives, Rd=numactives, Ru, D, numstates };

    enum obs { CASES, DEATHS, numobscolumns};

    enum excolumns {
#ifdef GOOGLE
      retail,	grocery,parks,	transit,	workplaces,	residential,
#endif
        REDUCTIONMEAN,
      R0,
                    numexcolumns};

    virtual unsigned k() const  { return numstates; }
    virtual unsigned n() const  { return numobscolumns; }

    std::string statelabel(unsigned i) const
    {
        static std::vector<std::string> r =
        { "E", "Ia", "Ip" ,"Is","Rd", "Ru","D" };
        return r[i];
    }

    enum params { v, ce, cs, ccoef, dcoef,
                  alpha, sigma, varsigma, gammas, mu, beta, tshift, ralpha,
#ifdef GOOGLE
                  cretail, firstgoogle=cretail,
                  cgrocery,cparks,	ctransit,	cworkplaces,cresidential,
                  lastgoogle = cresidential,
#endif
                  betafactor, mufactor, numparams};


    virtual dvector I(unsigned , const vector<double>& , const G& ) const
    {
        dvector ret(k());
        ret.setZero();
        return ret;
    }

    virtual dmatrix P(unsigned t, const vector<double>& params, const G& g ) const
    {
        dmatrix Pt(k(),k());
        Pt.setZero();

        Pt(E,Ia) = params[sigma] * params[alpha];
        Pt(E,Ip) = params[sigma] * (1-params[alpha]);
        Pt(Ip,Is) = params[varsigma];
        Pt(Ia,Ru) = params[gammas];
        Pt(Is,Rd) = params[gammas];
        double m = params[mu];
        if(t >= g.Ysize() - fineobs - g.estoffset)
            m *= params[mufactor];

        Pt(Is,D) = m;

        for(unsigned i=0; i<k(); i++)
        {
            double s = 0;
            for(unsigned j=0; j<k(); j++)
                s += Pt(i,j);
            Pt(i,i)=max(0.0,1-s);
        }
        return Pt.transpose();

    }

    virtual dmatrix hatB(unsigned t, const vector<double>& params, const G& g) const
    {
        dmatrix Bt(k(),k());
        Bt.setZero();

        double shiftedt = max(t-params[tshift],0.0);
        unsigned paqtl = static_cast<unsigned>(shiftedt);
        unsigned paqth = static_cast<unsigned>(shiftedt+1);
        double lweight = paqth - shiftedt;


        unsigned rgap = paqtl > 50 ? 50 : paqtl;

        unsigned rt = paqtl-rgap;
        double r0 = lweight*g.Z(rt,R0) + (1-lweight)*g.Z(rt+1,R0);

        for(;++rt <= paqtl;)
            r0 = params[ralpha] * lweight*g.Z(rt,R0) + (1-lweight)*g.Z(rt+1,R0)
                 + (1-params[ralpha]) * r0;

        double rl;
        double rr;
        double bet = params[beta];
#ifdef GOOGLE
        rl = bet;
        rr = bet;
        unsigned zp=0;
        for(unsigned i=firstgoogle; i<=lastgoogle; i++,zp++)
        {
            rl += params[i] * g.Z(paqtl,zp);
            rr += params[i] * g.Z(paqth,zp);
        }
        double b = (lweight*rl + (1-lweight)*rr)
              * r0;
#else
        rl = g.Z(paqtl,REDUCTIONMEAN);
        rr = g.Z(paqth,REDUCTIONMEAN);
        double b = (lweight*rl + (1-lweight)*rr)
              * r0 * bet;

#endif
        if(t >= g.Ysize() - fineobs - g.estoffset)
            b *= params[betafactor];


        Bt(Is,E) = b;
        Bt(Ip,E) = b;
        Bt(Ia,E) = b / 4;

        return Bt.transpose();
    }

    virtual double vb(const vector<double>& params ) const { return params[v]; }
    virtual double vp(const vector<double>& params, unsigned i) const
    {

       static int p[numstates ] =
      // E,  Ia,  Ip,  Is, Rd, Ru, D, numstates };
       { ce, ce,  ce,  cs, -1, -1, -1 };

        if(p[i]==-1)
            return numeric_limits<double>::infinity();
        else
            return params[p[i]];
    }


    virtual dvector gamma(unsigned /* t */, const vector<double>& params, const struct G& /*g*/) const
    {
        dvector ret(n());
        ret.setZero();
        ret[CASES] = params[ccoef];
        ret[DEATHS] = params[dcoef];
        return ret;
    }


    virtual dmatrix F(unsigned, const vector<double>& , const struct G& ) const
    {
       dmatrix ret(numobscolumns ,numstates);
       ret.setZero();

       ret(CASES, Is) = ret(CASES, D) = ret(CASES, Rd) =  1;
       ret(DEATHS,D) = 1;
       return ret;
    }

};


inline vector<double> weekadjust(const vector<double>& src, unsigned inercity, double threshold = 20, double outliertrsh = 30)
{

    if(src.size() < 8)
        throw "too small size to weekadjust";

    vector<double> dst;

    vector<double> lns;

    unsigned ignored = 0;
    for(unsigned i=0; i<src.size(); i++)
    {
        double x = max(i ? src[i]-src[i-1] : src[i],0.0);

        lns.push_back(x< numeric_limits<double>::epsilon()
                      ? numeric_limits<double>::quiet_NaN() : log(x));
        if(src[i] < threshold)
            ignored ++;
    }
    vector<double> cs(ignored,0);
    unsigned i=ignored;
    for(; i+3<lns.size(); i++)
    {
        if(lns[i]==numeric_limits<double>::quiet_NaN())
            cs.push_back(numeric_limits<double>::quiet_NaN());
        else
        {
            double s=0;
            unsigned n=0;
            for(unsigned j=i-3;j<=i+3; j++)
                if(!isnan(lns[j]))
                {
                    s += lns[j];
                    n++;
                }
            cs.push_back(n ? lns[i] - s / n : numeric_limits<double>::quiet_NaN());
        }
    }
    for(; i<lns.size(); i++)
        cs.push_back(cs[i-7]);

//    vector<double> acs;
    vector<double> ds;

    for(unsigned i=0; i<cs.size(); i++)
    {
        double s = 0;
        unsigned n =0;
        for(unsigned j=0; j<inercity; j++)
        {
            int k= i - 7 * j;
            if(k < 0)
                break;
            if(!isnan(cs[k]))
            {
                s+=cs[k];
                n++;
            }
        }
        double newc = n==0 ? 1 : exp(s/n);
//        acs.push_back(newc);

        ds.push_back((i ? src[i]-src[i-1] : src[i]) / newc);
    }

    double sumds = 0;
    for(int i=0; i<ds.size(); i++)
    {
        double s=0;
        double s2=0;
        unsigned n=0;
        constexpr int r = 3;
        for(int j=max(0,i-r);
            j<min(static_cast<int>(ds.size()),i+r); j++)
        {
            if(j!=i)
            {
                s += ds[j];
                s2 += ds[j]*ds[j];
                n++;
            }
        }
        double mean = s / n;
        double stdev = sqrt( s2 / n - mean * mean );
        if(ds[i] < 0 || (mean > outliertrsh && fabs(ds[i] - mean)/stdev > 5))
        {
            clog << "Outlier " << ds[i] << " at " << i << " replaced by " << mean << endl;
            ds[i] = mean;
        }
        sumds += ds[i];
    }

    double d = sumds / src[src.size()-1];

    double sumdst = 0;
    for(unsigned i=0; i<src.size(); i++)
    {
            sumdst += ds[i] / d;
            dst.push_back(sumdst);
//clog << src[i] << "," << dst[i] << endl;
    }
    return dst;
}

class ldatareader
{
public:

    time_t date2t(const string s)
    {
        struct tm ti;
        ti.tm_sec=0;
        ti.tm_min=0;
        ti.tm_hour=9;
        ti.tm_mday= stoi(s.substr(8,2));
        ti.tm_mon = stoi(s.substr(5,2))-1;
        ti.tm_year=stoi(s.substr(0,4))-1900;

        return mktime(&ti);
    }

    int date2int(const string s)
    {
        time_t res = date2t(s);
        int di = round((res + 86400 *0.5) / 86400);
        return di;
    }


    string t2date(time_t t)
    {
        struct tm * ptm = gmtime ( &t );
        ostringstream s;
        s << 1900 + ptm->tm_year << "-";
        if(ptm->tm_mon+1 < 10)
            s << "0";
        s << ptm->tm_mon+1 << "-";
        if(ptm->tm_mday < 10)
            s << "0";
        s << ptm->tm_mday;
        return s.str();
    }

//    ldatareader(const string& acountry) : country(acountry)
//    {}

    static constexpr double r0init = 4;
    static constexpr unsigned r0shift = 14;
    static constexpr unsigned forecastlen = 35;

    seirdata read(const string& googlefn, const string& country, const string& jhupath)
    {
        string firstdate = "2020-02-16";
        csv<','> cases(jhupath + "time_series_covid19_confirmed_global.csv" );
        csv<','> deaths(jhupath + "time_series_covid19_deaths_global.csv" );
        csv<','> mobility(googlefn);

        unsigned firstcol = 29;
        if(cases(0,firstcol) != "2/16/20")
            throw "JH changed format!";
        if(deaths(0,firstcol) != "2/16/20")
            throw "JH changed format!";
        if(mobility(1,0) != firstdate)
            throw "Vitek changed format!";

        seirdata res;
        res.ylabels.push_back("CASES");
        res.ylabels.push_back("DEATHS");


#ifdef GOOGLE
        res.zlabels.push_back("retail");
        res.zlabels.push_back("grocery");
        res.zlabels.push_back("parks");
        res.zlabels.push_back("transit");
        res.zlabels.push_back("workplaces");
        res.zlabels.push_back("residential");
#endif
        res.zlabels.push_back("REDUCTIONMEAN");
        res.zlabels.push_back("R0");

        unsigned j=1;
        for(; j<cases.r(); j++)
        {
            if(cases(j,1)==country && cases(j,0) == "")
                break;
        }

        if(j == cases.r())
            throw country + " not found!";

        if(cases(j,1) != deaths(j,1))
            throw "Country on different places";

        vector<double> c;
        vector<double> d;
        for(unsigned i=firstcol; i< cases.c(0); i++)
        {
            c.push_back(cases.getunsigned(j,i));
            d.push_back(deaths.getunsigned(j,i));
        }

//        vector<double> lastobs(numobscolumns,0);
//        res.z.resize(res.y.size());


        unsigned i=0;

        vector<double> r;
        bool any = false;
        for(; i<12; i++)
        {
            if(c[i] > 0)
                any = true;
            r.push_back(0);
        }

        for(; i< c.size(); i++)
        {
            if(c[i] > 0)
                any = true;
            if(!any)
                r.push_back(0);
            else
            {
                double num = c[i] - c[i-7];
                double den = c[i-5] - c[i-7-5];
                if(fabs(den) < 0.001)
                    r.push_back(r0init);
                else
                    r.push_back(num/den);
            }
        }

        vector<double> mr;

        constexpr unsigned mobcol = 7;
        double firstmob = mobility.getdouble(1,mobcol);
        int firstmobdate = date2int(mobility(1,0));

        for(unsigned i=1; i<mobility.r()-1; i++)
        {
            for(unsigned j=0; j<7; j++)
                mr.push_back((mobility.getdouble(i,mobcol) * (1 - j / 7.0)
                             +mobility.getdouble(i+1,mobcol) * j / 7.0)/firstmob);
            if(date2int(mobility(i,0))-firstmobdate != 7 * (i-1))
                throw "mobility data not ascening";
        }
#ifdef GOOGLE

        constexpr unsigned ngoogles = 6;
        vector<vector<double>> ggs(ngoogles);

        for(unsigned i=1; i<mobility.r()-1; i++)
        {
            for(unsigned j=0; j<7; j++)
                for(unsigned k=0; k<ngoogles; k++)
                    ggs[k].push_back((mobility.getdouble(i,1+k) * (1 - j / 7.0)
                             +mobility.getdouble(i+1,1+k) * j / 7.0));
        }
#endif
        vector<double> r0;
        for(unsigned i=0; i+r0shift < r.size(); i++)
            r0.push_back(r[i+r0shift] / mr[i]);

        vector<double> cs = weekadjust(c,4,100);

        for(unsigned i=0; i<cs.size(); i++)
            res.y.push_back(dv({cs[i],d[i]}));

        for(unsigned i=0; i<res.y.size()+forecastlen; i++)
            res.z.push_back(dv({
#ifdef GOOGLE
                             ggs[0][min(i,static_cast<unsigned>(ggs.size()-1))],
                               ggs[1][min(i,static_cast<unsigned>(ggs.size()-1))],
                               ggs[2][min(i,static_cast<unsigned>(ggs.size()-1))],
                               ggs[3][min(i,static_cast<unsigned>(ggs.size()-1))],
                               ggs[4][min(i,static_cast<unsigned>(ggs.size()-1))],
                               ggs[5][min(i,static_cast<unsigned>(ggs.size()-1))],
#endif
                             mr[min(i,static_cast<unsigned>(mr.size()-1))],
                             r0[min(i,static_cast<unsigned>(r0.size()-1))]}));
        time_t t = date2t(firstdate);
        for(unsigned i=0; i<res.y.size()+forecastlen; i++, t+= 24*3600)
            res.dates.push_back(t2date(t));
        res.lag = 0;
        return res;
    }
};


template <estimationmethod M=emwls>
class estlseir : public lseir, public estimableseirfilter<M>
{
public:
    virtual uncertain X0(const vector<double>& params) const
    {
        return fx0;
    }

    estlseir( uncertain x0 ) : fx0(x0) {}
    uncertain fx0;
};


#endif // LSEIR_HPP
