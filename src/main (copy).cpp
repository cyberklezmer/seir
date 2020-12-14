#include <vector>
#include "emodel.hpp"
#include "orpp.hpp"
#include "seir.hpp"
#include <cmath>

using namespace std;
using namespace orpp;

struct paramval
{
    string n;
    double v;
};

//#define POSITIVITY
#define EXPAPPROX
#define DONOTSUM
//#define QUARANTINERATIO
//#define VACATIONDUMMY
#define TEMPER
//#define SHIFT
//#define FACTORPAQ
//#define MU
//#define FAMILYCONTACTS
//#define BETAADD
#define ARATE
#define ASYMP
//#define CONSTIOTA
//#define PIECEWISEIOTA
//#define SETA
//#define PIECEWISEETA
#define GAMMAS

#define PIECEWISETHETA
//#define SCURVE

//#define K
//#define KACTIVE
//#define K2

static constexpr unsigned horizon=246;
static constexpr unsigned dusekpenalty=00;

vector<unsigned> thetafrns({20,60,130,160,180,220});

struct paraminit { string n;	double v;	bool filtered; };

vector<paraminit> initvals = {
    {"marchimportrate",9.97028,true},
    {"beta",0.659874,true},
    {"phi",1.21973e-20,false},
    {"asymprate",0.545619,true},
    {"omega1",1.10711,false},
    {"omega2",0.782054,false},
    {"rho0",0.0106261,true},
    {"rho1",0.00453292,true},
    {"rho2",0.0195615,true},
    {"rho3",0.0272741,true},
    {"rho4",0.0273509,true},
    {"rho5",0.010006,true},
    {"rho6",0.0193714,true},
    {"eta",0.682174,false},
    {"gammaa",77,true},
    {"gammas",60,true},
    {"gammad",0,true},
};


enum eparams {
    emarchimportrate,
               ebeta,
#ifdef BETAADD
    ebetaadd,
#endif

#ifdef CONSTIOTA
    eiotaconst,
#endif
#ifdef TEMPER
    ephi,
#endif

#ifdef VACATIONDUMMY
    evdummy,
#endif

#ifdef SHIFT
    eshift,
#endif
#ifdef ARATE
    easymprate,
#endif

#ifdef MU
    emu,
#endif
#ifdef QUARANTINERATIO
    eqratio,
#endif

#ifdef PIECEWISEIOTA
    eiota0,eiota1,eiota2,eiota3,eiota4,eiota5,
    eiota6, eiota7,eiota8,eiota9,eiota10, elastiota = eiota10,
#else
#ifdef FACTORPAQ
    eomega0, eomegas,eomegag,
#else
#ifdef FAMILYCONTACTS
    efcontacts,
#endif
    eomega1,  eomega2,
#endif
#endif // PIECEWISEIOTA

#ifdef PIECEWISETHETA
               erho0,erho1, erho2, erho3, erho4, erho5,
               erho6, elastrho = erho6,
#else
#ifdef SCURVE
        erhoinit, erhosize, erhomid, erhok,
#else
        erho,
#endif
#endif
#ifdef POSITIVITY
    ethetaposcoef,eetaposcoef,
#endif

#ifdef PIECEWISEETA
               eeta0,eeta1,eeta2,eeta3,eeta4,eeta5,
               eeta6, elasteta = eeta6,
#else
#ifdef SETA
                  eetamin,
        eetasize, eetamid, eetak,
#else
               eeta,
#endif
#endif

#ifdef K
    ek,
#endif
#ifdef KACTIVE
    ekactive,
#endif

#ifdef K2
    ek2,
#endif

#ifdef GAMMAS
               egammaa,egammas,egammad,
#endif

      enumparams};


enum eseirstates {eseire,
                  eseiria, eseiris, eseirin, eseirjs, eseirjn, eseirr, eseird,
                  eseirqe,

                  eseiriada, eseirinda,
                  eseirisda, eseirjsda, eseirjnda, eseirrda, eseirdda,

                  eseirisds, eseirjsds,eseirrds, eseirdds,

                  eseirnumstates
                 };

class pwfn
{
    vector<unsigned> frs;
    vector<unsigned> trsf;
public:
    pwfn(const vector<unsigned>& frontiers, int numpars=-1)
        : frs(frontiers),
          trsf(frs.size() == 0 ? 1 : frontiers[frontiers.size()-1]+1)
    {
        assert(numpars == -1 || frontiers.size()+1 == numpars);
        unsigned j=0;
        unsigned i=0;
        for(; i<frontiers.size(); i++)
        {
            assert(frontiers[i] > j);
            for(;j<frontiers[i];j++)
                trsf[j]=i;
        }
        trsf[j] = i;
    }
    double operator () (unsigned t, const vector<double>& p, unsigned offset=0)
    {
        unsigned numfrs = trsf[trsf.size()-1];
        assert(offset + numfrs < p.size());

        if(t>=trsf.size())
            return p[offset+numfrs];
        unsigned i1 = trsf[t];
        double t1 = i1 ? frs[i1-1] : 0;
        double t2 = frs[i1];
        double w = (t2 - t) / (t2-t1);
        return p[offset + i1] * w + p[offset + i1+1] * (1-w);
    }
};



#ifdef PIECEWISETHETA
    pwfn thetafn(thetafrns, elastrho-erho0+1);
#endif

#ifdef PIECEWISEETA
    pwfn etafn(thetafrns, elasteta-eeta0+1);
#endif


inline string seirstatelabel(unsigned i)
{
    static const string l[] = {
                   "E",
                   "Ia", "Is", "In", "Js", "Jn", "R", "D",
                   "QE",

                   "Iada", "Inda",
                   "Isda", "Jsda","Jnda", "Rda", "Dda",

                   "Isds", "Jsds","Rds", "Dds"
    };
    return l[i];
}

inline double zofunction(double x, double c = 0.05)
{
    if(x < c)
        return c* exp((x-c) / c);
    else if(x > 1-c)
        return 1 - c* exp((1-x-c) / c);
    else
        return x;
}

class czseir : public seirmodel
{
public:
    enum eobs {
#ifdef ASYMP
        eRA, eRS,
#else
        eR,
#endif
        eD,
#ifdef HOSPITALS
        eH, eHC,
#endif
        enumobs};
    enum eexp { eI, econtacts, egammared, esred, etests, edayadjust,
                enumexp};
private:
    static matrix getH()
    {
        matrix ret(czseir::enumobs,eseirnumstates,0);
#ifdef ASYMP       
        for(unsigned i=eseiriada; i<=eseirdda; i++)
            ret(eRA,i) = 1;
        for(unsigned i=eseirisds; i<=eseirdds; i++)
            ret(eRS,i) = 1;
        ret(eRS,eseird)=1;
#else
        for(unsigned i=eseiriada; i<eseirnumstates; i++)
            ret(eR,i) = 1;
        ret(eR,eseird)=1;
#endif
        ret(eD,eseird)=1;
        ret(eD,eseirdda)=1;
        ret(eD,eseirdds)=1;
        return ret;
    }
public:
    czseir() : seirmodel(getH()) {}

    matrix P(unsigned  t, const vector<double>& params, const G& g) const
    {
        double asymprate;
        expparams p;
#ifdef ARATE
        asymprate = params[easymprate];
#else
        asymprate = p.asymptomatic_rate;
#endif
        double sigma = p.sigma; //ksigma];
        vector<vector<double>> ret(eseirnumstates,vector<double>(eseirnumstates,0));
        ret[eseire][eseirin] = sigma * asymprate;
        ret[eseire][eseiria] = sigma * (1-asymprate);

        double thetas = getthetas(t,params,g) * g.z[t][czseir::edayadjust];
        double thetaa = getthetaa(t,params,g) * g.z[t][czseir::edayadjust];
        ret[eseire][eseirqe] = thetaa;
        ret[eseirqe][eseirinda] = sigma * asymprate;
        ret[eseirqe][eseiriada] = sigma * (1-asymprate);

        ret[eseirin][eseirjn] = p.delta_n;
        ret[eseirjn][eseirr] = p.gamma_In;
        ret[eseiria][eseiris] =
                p.symptoms_manifest_rate;

        ret[eseiris][eseirjs] = p.delta_s;
        ret[eseirjs][eseirr] = p.gamma_Is;
#ifdef MU
        ret[eseirjs][eseird] = params[emu];
#else
        ret[eseirjs][eseird] = p.mu;
#endif
        ret[eseiria][eseiriada] = thetaa;
        ret[eseirin][eseirinda] = thetaa;
        ret[eseiris][eseirisds] = thetas,

        ret[eseirinda][eseirjnda] = p.delta_n;
        ret[eseirjnda][eseirrda] = p.gamma_In;
        ret[eseiriada][eseirisda] = p.symptoms_manifest_rate;
        ret[eseirisda][eseirjsda] = p.delta_s;
        ret[eseirisds][eseirjsds] = p.delta_s;
        ret[eseirjsda][eseirrda] = p.gamma_Is;
        ret[eseirjsds][eseirrds] = p.gamma_Is;
#ifdef MU
        ret[eseirjsda][eseirdda] = params[emu];
        ret[eseirjsds][eseirdds] = params[emu];
#else
        ret[eseirjsda][eseirdda] = p.mu;
        ret[eseirjsds][eseirdds] = p.mu;
#endif
        for(unsigned i=0; i<eseirnumstates; i++)
        {
            double s=0;
            for(unsigned j=0; j<eseirnumstates; j++)
                if(i != j)
                    s+=ret[i][j];
            if(s>1)
            {
                if(s > 1.05)
                    throw "s > 1!";
                ret[i][i] = 0;
            }
            else
                ret[i][i] = 1.0-s;
        }
        return ret;
    }
    vector<double> I(unsigned  t, const vector<double>& params, const G& g) const
    {
        vector<double> ret(eseirnumstates,0);
        double coef = t < 30 ? params[emarchimportrate] : 1.0 / (1.0 - 0.179);

        ret[eseire]= coef * g.z[t][eI];
        return ret;
    }


private:
    virtual vector<double> getiota(unsigned t, const vector<double>& p, const G& g) const = 0;
    virtual double getthetas(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual double getthetaa(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual string statelabel(unsigned i) const { return seirstatelabel(i); }
public:
    double numantibodies(const evalresult& r, unsigned t) const
    {
        double s=0;
        vector<double> est = r.est[t].x();
        for(unsigned i=0; i< eseirnumstates; i++)
        {
            switch(i)
            {
//                case eseire:
//                case eseirqe:
                case eseird:
                case eseirdda:
                case eseirdds:
                    break;
                default:
                    s += est[i];
            }
        }
        return s;
    }

    double numactive(const G& g, unsigned t) const
    {
        double s=0;
        vector<double> est = g.est[t];
        for(unsigned i=0; i<= eseirin ; i++)
            s += est[i];
        return s;
    }


    virtual double contrastaddition(const evalresult& r) const
    {
        unsigned t = 70;

        double s = this->numantibodies(r,t);
        // 107 of 26500, population 10693000
        double var = 4173.0 * 4173.0;
        double mu = 43175;
        double ll =
                -0.5 * ::log(2*3.1415926)
                    -0.5 * ::log(var)
                    -0.5 * (s-mu)*(s-mu) / var;
        return dusekpenalty * ll;
    }


    matrix J(unsigned t, const vector<double>& params, const G& g) const
    {
        vector<double> iot = getiota(t,params,g);

        matrix ret(k(),k(),0);
        for(unsigned i=0; i<eseirnumstates; i++)
        {
            assert(iot[i] >= 0);
            ret(i,0) = iot[i];
        }
        return ret;
    };
    vector<double> X0(const vector<double>& , const G&) const
    {
        vector<double> ret(k(),0);
//        ret[eseire] = params[ealpha];
        return ret;
    }
    virtual matrix Gamma(unsigned /* t */, const vector<double>& params, const G& ) const
    {
        double gammad;
        double gammaa;
        double gammas;
#ifdef GAMMAS
        gammaa = params[egammaa];
        gammas = params[egammas];
        gammad = params[egammad];
#else
        gammaa = 0;
        gammas = 0;
        gammad = 0;
#endif

        matrix ret = matrix(m(),k(),0);
#ifdef ASYMP
        ret(eRS,eseiris) = gammas;
        ret(eRA,eseire)=ret(eRA,eseiria)=
                ret(eRA,eseirqe)=ret(eRA,eseirin)=gammaa;
#else
        ret(eR,eseiris) = gammas;
#endif
        ret(eD,eseirjsda) = gammad;
        ret(eD,eseirjsds) = gammad;
        ret(eD,eseirjs) = gammad;
        return ret;
    }
};


class obsseir : public czseir
{
public:

    vector<double> getiota(unsigned at, const vector<double>& p, const G& g) const
    {
        unsigned t1, t2;
        double w;
#ifdef SHIFT
        double tt = at-p[eshift];
        t1 = max(0.0,floor(tt));
        t2 = min(t1+1,at);
        w = t2-tt;
#else
        t1=t2=at;
        w=1;
#endif


        assert(p.size()==enumparams);
        double ba=0;
#ifdef BETAADD
        ba = p[ebetaadd];
#endif
#ifdef PIECEWISEIOTA
        double preiota;
        unsigned numiotas = elastiota + 1 - eiota0;
        unsigned period = horizon / numiotas;
        unsigned i = min(at / period, numiotas - 1);
        if(i == numiotas - 1)
            preiota = p[eiota0 + i];
        else
        {
            double w = (at % period) / static_cast<double>(period);
            preiota = (1-w) * p[eiota0+i] + w * p[eiota0+i+1];
        }
#else
#ifdef FACTORPAQ
double preiota = // exp(-p[eomega0] - p[eomegag]*g.z[t1][egammared] - p[eomegas]*g.z[t1][esred]);
                zofunction(1.0-p[eomegas]*(w* g.z[t1][esred] + (1-w)* g.z[t2][esred]))
                * zofunction(p[eomega0] -p[eomegag]*(w* g.z[t1][egammared] + (1-w)* g.z[t2][egammared]));
#else
        double c = w*g.z[t1][econtacts]+(1-w)*g.z[t2][econtacts];
#ifdef FAMILYCONTACTS
        c -= p[efcontacts];
#endif
#ifdef VACATIONDUMMY
        if(t1>128 && t1 < 189)
            c *= p[evdummy];
#endif

        double gammared = w*g.z[t1][egammared]+(1-w)*g.z[t2][egammared];
#ifdef EXPAPPROX
        double preiota =c * exp(-p[eomega1]*gammared - p[eomega2]*gammared*gammared);
#else
        double preiota =c * pow(1-p[eomega1]*gammared, p[eomega2]);
#endif
        if(preiota  > 1.5)
            throw("preiota");

#endif
#endif // PIECEWISEIOTA
#ifdef CONSTIOTA
        preiota += p[eiotaconst];
#endif
        double d;
#ifdef TEMPER
        const unsigned t0 = -24;
        d = (1 - p[ephi] * cos(-t0 / 365.0 * 2* 3.141592)) + p[ephi] * cos((at-t0) / 365.0 * 2 * 3.141592 );
#else
        d = 1;
#endif
        assert(preiota >= 0);
        vector<double> iota(eseirnumstates,0);
        iota[eseiria] = iota[eseiris] = (p[ebeta]+ba)*d * preiota;
        iota[eseirin] = (p[ebeta]+ba)* d *preiota;
        assert(iota[eseiria] >= 0);
        assert(iota[eseirin] >= 0);
        assert(iota[eseiris] >= 0);
        return iota;
    }

    double weekave(unsigned t, const G& g) const
    {
        unsigned former = t>7 ? t-7 : 0;
        double numrep;
        if(t == 0)
            numrep=0;
        else
        {
#ifdef ASYMP
            numrep = (g.y[t][eRA]+g.y[t][eRS]-g.y[former][eRA]-g.y[former][eRS])/ (t-former);
#else
            numrep = (g.y[t][eR]-g.y[former][eR])/ (t-former);
#endif
        }
        return numrep;
    }

    virtual double getthetas(unsigned t, const vector<double>& params, const G& g) const
    {
        double rh;
#ifdef PIECEWISEETA
/*
        unsigned numetas = elasteta + 1 - eeta0;
        unsigned period = horizon / numetas;
        unsigned i = min(t / period, numetas - 1);
        if(i == numetas - 1)
            rh = params[eeta0 + i];
        else
        {
            double w = (t % period) / static_cast<double>(period);
            rh = (1-w) * params[eeta0+i] + w * params[eeta0+i+1];
        }
*/
        rh = etafn(t,params, eeta0);
#else // PIECEWISEETA
#ifdef SETA
        double eta = params[eetamin]
                + params[eetasize] / (1.0 + exp(-params[eetak] * (t -params[eetamid] )));
#else
        double eta = params[eeta];
#endif
#ifdef DONOTSUM
        rh = eta;
#else
        rh = getthetaa(t,params,g) + eta;
#endif
#endif
#ifdef POSITIVITY
      double positivity = g.z[t][etests] ? weekave(t,g) / g.z[t][etests] : 0;
      rh *= 1-params[eetaposcoef] * positivity;
#endif

        return  zofunction(rh / 0.7,0.000001) *0.7;
    }


    virtual double getthetaa(unsigned t, const vector<double>& params, const G& g) const
    {
        double qratio;
   #ifdef QUARANTINERATIO
        qratio = params[eqratio];
   #else
        qratio = 1;
   #endif

        double rh;
#ifdef PIECEWISETHETA
/*        unsigned numrhos = elastrho + 1 - erho0;
        unsigned period = horizon / numrhos;
        unsigned i = min(t / period ,numrhos - 1);
        if(i == numrhos - 1)
            rh = params[erho0 + i];
        else
        {
            double w = (t % period) / static_cast<double>(period);
            rh = (1-w) * params[erho0+i] + w * params[erho0+i+1];
        }*/
        rh = thetafn(t,params, erho0);
#else
#ifdef SCURVE
       rh = params[erhoinit]+params[erhosize] / (1.0 + exp(-params[erhok] * (t -params[erhomid] )));
#else
       rh = params[erho];
#endif
#endif
       double k = 1e10;
#ifdef K
       if(t < 100)
          k = params[ek];
#endif
#ifdef K2
       if(t >= 100)
          k = params[ek2];

#endif
       double c;
#ifdef KACTIVE
     double na= numactive(g,t);
     c = na < params[ekactive] ? rh : params[ekactive] / na * rh
                                 + (na - params[ekactive]) / na * rh * qratio;
#else
      double numrep = weekave(t,g);

      c = numrep < k ? rh : ((1-qratio) * rh + qratio * rh * k / numrep);
#endif
#ifdef POSITIVITY
      double positivity = g.z[t][etests] ? weekave(t,g) / g.z[t][etests] : 0;
      c *= 1-params[ethetaposcoef] * positivity;
#endif
        return zofunction(c/0.7,0.000001) * 0.7; // max(0.0, min(0.7, c));
    }
};

void fillsi(seirmodel::data& si, csv<','>&c, unsigned lag=0)
{
    bool rfin = false;
    vector<unsigned> trsf(czseir::enumobs);
#ifdef ASYMP
    trsf[czseir::eRA] = 0;
    trsf[czseir::eRS] = 1;
#else
    trsf[czseir::eR] = 2;
#endif
    trsf[czseir::eD] = 3;
#ifdef HOSPITALS
    trsf[czseir::eH] = 4;
    trsf[czseir::eHC] = 5;
#endif
    const unsigned firstexp = 6;
    for(unsigned j=0; j<czseir::enumobs; j++)
    {
        si.ylabels.push_back(c(0,1+trsf[j]));
    }

    for(unsigned j=0; j<czseir::enumexp; j++)
        si.zlabels.push_back(c(0,firstexp+j));

    for(unsigned i=lag; i< c.r()-1; i++)
    {
        si.dates.push_back(c(i+1,0));
        if(c(i+1,1) != "" && !rfin)
        {
            vector<double> obs;
            for(unsigned j=0; j<czseir::enumobs; j++)
                obs.push_back(c.getunsigned(i+1,1+trsf[j]));
            si.y.push_back(obs);
        }
        else
            rfin = true;
        vector<double> exp;
        for(unsigned j=0; j<czseir::enumexp; j++)
            exp.push_back(c.getdouble(i+1,1+firstexp+j));
        si.z.push_back(exp);
    }
}



void seir()
{
    csv<','> c("../input/cz.csv");
    seirmodel::data si;

    fillsi(si,c,0);
    assert(horizon <= si.y.size());
    si.y.resize(horizon);
    si.z.resize(min(100UL+horizon,si.z.size()));

    auto siest = si;





    vector<paraminfo> params
            = {
                 paraminfo("marchimportrate",3,0,10),
                 paraminfo("beta",2.08,0,8),
#ifdef BETAADD
              paraminfo("betaadd",0,0,3),
#endif
#ifdef CONSTIOTA
        paraminfo("constiota",0,0,0.2),
#endif
    #ifdef TEMPER
        paraminfo("phi",0,0,0.45),
    #endif
    #ifdef VACATIONDUMMY
        paraminfo("vdummy",1,0,1),
    #endif

    #ifdef SHIFT
            paraminfo("shift",5,0,7),
    #endif

#ifdef ARATE
    paraminfo("asymprate", 0.4, 0.15, 0.55),
#endif

#ifdef MU
    paraminfo("mu", 0.00131326350277332, 0, 0.02),
#endif

#ifdef QUARANTINERATIO
        paraminfo("qratio", 0.1, 0, 1),
#endif

#ifdef PIECEWISEIOTA
        paraminfo("iota0",3.28,0,5),
        paraminfo("iota1",0.2,0,5),
        paraminfo("iota2",0.3,0,5),
        paraminfo("iota3",0.3,0,5),
        paraminfo("iota4",0.4,0,5),
        paraminfo("iota5",0.4,0,5),
        paraminfo("iota6",0.5,0,5),
        paraminfo("iota7",0.5,0,5),
        paraminfo("iota8",0.6,0,5),
        paraminfo("iota9",0.6,0,5),
        paraminfo("iota10",0.7,0,5),
#else
#ifdef FACTORPAQ
        paraminfo("omega0", 1, 0, 0.5),
        paraminfo("omegas", 1, 0, 1),
        paraminfo("omegag", 1, 0, 1),
#else
    #ifdef FAMILYCONTACTS
                    paraminfo("fcontacts",0.0,0,0.3),
    #endif
                paraminfo("omega1", 2, 0, 10),
                paraminfo("omega2", 0, 0, 10),

#endif
#endif // PIECEWISEIOTA

#ifdef PIECEWISETHETA
                paraminfo("rho0", 0, 0, 0.7),
                paraminfo("rho1", 0.1, 0, 0.7),
                paraminfo("rho2", 0.1, 0, 0.7),
                paraminfo("rho3", 0.2, 0, 0.7),
                paraminfo("rho4", 0.2, 0, 0.7),
                paraminfo("rho5", 0.3, 0, 0.7),
                paraminfo("rho6", 0.3, 0, 0.7),
#else
#ifdef SCURVE
            paraminfo("rhoinit", 0, 0, 0.3),
            paraminfo("rhosize", 0.10, 0, 0.7),
            paraminfo("rhomid", 80, 0, 300),
            paraminfo("rhok", 0.03, 0, 1),
#else
            paraminfo("rho", 0.21, 0, 0.7),
#endif
#endif

#ifdef POSITIVITY
        paraminfo("thetaposcoef", 1, 0, 2),
        paraminfo("etaposcoef", 1, 0, 2),
#endif

#ifdef PIECEWISEETA
                paraminfo("eta0", 0.032575, 0, 0.7),
                paraminfo("eta1", 0.032575, 0,0.7),
                paraminfo("eta2", 0.032575, 0, 0.7),
                paraminfo("eta3", 0.032575, 0, 0.7),
                paraminfo("eta4", 0.032575, 0, 0.7),
                paraminfo("eta5", 0.032575, 0, 0.7),
                paraminfo("eta6", 0.032575, 0, 0.7),
#else
#ifdef SETA
        paraminfo("etamin", 0.1, 0, 0.5),
        paraminfo("etasize", 0.1, 0, 0.3),
        paraminfo("etamid", 100, 0, 300),
        paraminfo("etak", 0.06, 0, 1),
#else
        paraminfo("eta", 0.3, 0, 0.7),
#endif
#endif

#if defined(K)
             paraminfo("k", 100, 0, 260),
#else
#ifdef KACTIVE
                    paraminfo("kactive", 5000, 0, 100000),
#endif
#endif
#ifdef K2
               paraminfo("k2", 500, 0, 2000),
#endif

#ifdef GAMMAS
                paraminfo("gammaa",0, 0 , 100 ),
                paraminfo("gammas",0, 0 , 100 ),
                paraminfo("gammad",0,0, 100),
#endif

              };

assert(params.size()==enumparams);

    vector<bool> filter(params.size(),false);
//    filter[eomega0]=filter[eomega1]=filter[efcontacts]=true;
//    for(unsigned i=erho0; i<=elastrho; i++)
//        filter[i] = true;
//    for(unsigned i=erho0; i<=elastrho; i++)
//        filter[i] = true;

    obsseir s;

    vector<double> rp(enumparams);

    for(unsigned i=0; i<params.size(); i++)
    {
        rp[i] = params[i].initial;
        for(unsigned j=0; j<initvals.size(); j++)
            if(params[i].name == initvals[j].n)
            {
                rp[i]=params[i].initial = initvals[j].v;
                filter[i] = initvals[j].filtered;
                cout << params[i].name << " set filter to " << filter[i] << endl;
            }
    }

    if(1)
    {
        if(0)
        {
            uncertain res;
            cout << "ll= " << s.estimate(params,siest,res,filter) << endl;
            for(unsigned i=0; i<params.size(); i++)
            {
                cout << params[i].name << " ";
                cout << res.x()[i] << "(" << sqrt(res.var()(i,i)) << ")" << endl;
                params[i].initial = res.x()[i];
            }
            rp = res.x();
            for(unsigned i=0; i<params.size(); i++)
            {
                cout << "{\"" << params[i].name << "\","
                    << rp[i] << "," << (filter[i] ? "true" : "false") << "},";
                cout << endl;
            }
        }

        csvout<double,','>(clog, rp);
        czseir::evalresult r = s.eval<true>(rp,si, si.z.size()-si.y.size());
        ofstream o("output.csv");
        r.output(o,si);
        clog << "ll/n=" << r.contrast / r.contrasts.size() << endl
             << "additional contrast = " << s.contrastaddition(r) << endl
             << "Duseks survey = " << s.numantibodies(r,70) << endl
             << endl << endl;
        if(1) // confidence of reported
        {
            cout << "stddevs R" << endl;
            for(unsigned i=0; i<r.pred.size(); i++)
            {
                matrix v = r.pred[i].var();
                cout << sqrt(v(eseirnumstates+czseir::eRA,eseirnumstates+czseir::eRA)
                        + v(eseirnumstates+czseir::eRS,eseirnumstates+czseir::eRS)
                        + 2 * v(eseirnumstates+czseir::eRA,eseirnumstates+czseir::eRS)) << endl; ;
            }
            throw;
        }

        if(0) // R
        {
            clog << "rho, R,lastval" << endl;
             for(unsigned t=0; t<si.y.size(); t++)
            {
                double R = s.R(t,r);
                double Rcp = s.Rcp(t,r,4);
                clog << s.rho(t,r,eseirin+1) << "," << R << "," << Rcp << endl;
            }
            throw;
        }

        if(0) // calculator
        {
            matrix T(eseirin+1,eseirin+1,0);
            for(unsigned i=0; i<7; i++)
                T = T + r.Ts[si.y.size()-1-i].submatrix(0,0,eseirin+1,eseirin+1);
            T = (1.0 / 7.0) * T;
            double gamma = T(0,1);
            double c=si.z[si.z.size()-1][czseir::econtacts];
            cout << "C,rho" << endl;
            for(unsigned i=30; i<100; i++)
            {
                double factor = (i/100.0) / c;
                T(0,1)=T(0,2)=T(0,3)=gamma*factor;
                cout << i/100.0 << "," << range(T) << endl;
            }
            throw;
        }



        clog <<
#ifdef ASYMP
                "PA,predPA, predPAtrend, resA, resAtrend,"
                "PS,predPS, predPStrend, resS, resStrend,"
#else
                "P,predP, predPtrend, res, restrend,,,,,,"
#endif
                "D,predD, predDtrend, , , iota, eta, theteq,rho,rtrend,"
                "R,  ,contrast" << endl;
        for(unsigned i=0; i<si.y.size(); i++)
        {
            double D = si.y[i][czseir::eD];
#ifdef ASYMP
            double PA = si.y[i][czseir::eRA];
            double PS = si.y[i][czseir::eRS];
            double predPAtrend;
            double predPStrend;
#else
            double P = si.y[i][czseir::eR];
            double predPtrend;
#endif
            double predDtrend;
            matrix trendm;
            if(i>15)
            {
                unsigned ds = min(i-15,7U);
                auto six = si;
                six.y.resize(i-ds+1);
                czseir::evalresult rx =
                        s.eval<true>(rp,six,ds);
#ifdef ASYMP
                predPAtrend = rx.pred[i].x()[eseirnumstates+czseir::eRA];
                predPStrend = rx.pred[i].x()[eseirnumstates+czseir::eRS];
#else
                predPtrend = rx.pred[i].x()[eseirnumstates+czseir::eR];
#endif
                predDtrend = rx.pred[i].x()[eseirnumstates+czseir::eD];
                trendm = rx.Ts[i].submatrix(0,0,eseirin+1,eseirin+1);
            }
            else
            {
#ifdef ASYMP
                predPAtrend = PA;
                predPStrend = PS;
#else
                predPtrend = P;
#endif
                predDtrend = D;
                trendm = r.Ts[i].submatrix(0,0,eseirin+1,eseirin+1);
            }
#ifdef ASYMP
            double predPA = r.pred[i].x()[eseirnumstates+czseir::eRA];
            double predPS = r.pred[i].x()[eseirnumstates+czseir::eRS];
#else
            double predP = r.pred[i].x()[eseirnumstates+czseir::eR];
#endif
            double predD = r.pred[i].x()[eseirnumstates+czseir::eD];
//            double predPtrend = r2.pred[i].x()[eseirnumstates+czseir::eT];

            double rng = range(r.Ts[i].submatrix(0,0,eseirin+1,eseirin+1));

            double rn = s.R(i,r);
            double rcp = s.Rcp(i,r,4);

            clog
#ifdef ASYMP
                    << PA << "," << predPA << "," << predPAtrend << ","
                    << (predPA - PA)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eRA,eseirnumstates+czseir::eRA)) << ","
                    << predPAtrend - PA << ","
                    << PS << "," << predPS << "," << predPStrend << ","
                    << (predPS - PS)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eRS,eseirnumstates+czseir::eRS)) << ","
                    << predPStrend - PS << ","
#else
                 << P << "," << predP << "," << predPtrend << ","
                 << (predP - P)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eR,eseirnumstates+czseir::eR)) << ","
                 << predPtrend - P << ",,,,,,"
#endif
                 << D << "," << predD << "," << predDtrend << ","
                         << predD - D << ","
                         << predDtrend - D << ","
                 << r.Ts[i](eseire,eseiris) << ","
                 << r.Ts[i](eseirisds,eseiris) << ","
                  << r.Ts[i](eseirqe,eseire) << ","
                 << rng << ","
                << range(trendm) << ","
                << rn << "," << rcp;
            if(i<r.contrasts.size())
                clog << "," << r.contrasts[i];
            clog << endl;
        }
        for(unsigned i=0; i<params.size(); i++)
        {
            cout << "{\"" << params[i].name << "\","
                << rp[i] << "," << (filter[i] ? "true" : "false") << "},";
//            cout << "(" << sqrt(r.var()(i,i)) << ")"
            cout << endl;
        }


    }
    else
    {
/*        vector<double> res;
        filter[egamma]=true;
        cout << "ll= " << s.estimatetrend(params,res,filter) << endl;
        for(unsigned i=0; i<params.size(); i++)
        {
            cout << params[i].name << " ";
            cout << res[i] << endl;
        }
        vector<vector<double>> r = s.evaltrend(res).second;
        for(unsigned i=0; i<r.size(); i++)
        {
            cout << r[i][0] << "," << si.obs[i][0] << endl;
        }*/
    }
}


int main()
{
    double startt = time(0);

    try
    {
        seir();
    }
    catch (std::exception& e) {
        std::cerr << e.what() << endl;
        return 1;
    }
    catch (const char* m)
    {
           cerr << m << endl;
           return 1;
    }
    catch (const string& m)
    {
           cerr << m << endl;
           return 1;
    }
    catch(...)
    {
        cerr << "unknown exception" << endl;
        return 1;
    }


    clog << "time: " << time(0)-startt << " seconds." << endl;

    return 0;
}
