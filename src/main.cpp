#include "orpp.hpp"
#include "orpp/matrices.hpp"
#include "orpp/csv.hpp"




#include <vector>
#include "emodel.hpp"
#include "seir.hpp"
#include <cmath>

using namespace std;
using namespace orpp;

struct paramval
{
    string n;
    double v;
};

struct paraminit { string n;	double v;	bool filtered; };


//#define GENERALTEST

// #define TRV1
// #define TRV2
//#define NORMAL
#define WORKING
#ifdef NORMAL



static constexpr seirmodel::contrasttype ct = seirmodel::ectwls ; //ectll // ectwls; // ;:;
static constexpr bool nopred = false;
static constexpr bool estimate = true;
static constexpr bool sensanal = true;

static constexpr double timest = 20*60;

static constexpr unsigned horizon=295;
static constexpr unsigned dusekpenalty=0;

#define NOEXTRAMIMPORT
#define DEADFROMIS
#define PDETTHETA
#define PDET
#define PDETCOEF
#define PDETONLYETA
//#define PARKS
#define IMUNITY
#define DAKTELA
#define VARFACTOR
//#define POSITIVITY
#define EXPAPPROX
//#define LINAPPROX
//#define DONOTSUM
//#define QUARANTINERATIO
//#define VACATIONDUMMY
//#define TEMPER
#define SHIFT
//#define FACTORPAQ
#define MU
//#define FAMILYCONTACTS
//#define BETAADD
//#define ARATE
#define ASYMP
//#define CONSTIOTA
//#define PIECEWISEIOTA
//#define SETA
//#define PIECEWISEETA
#define GAMMAS
//#define PIECEWISETHETA
//#define THREEPIECETHETA
//#define SCURVE

//#define K
//#define KACTIVE
//#define K2
#define WEIGHTVACCATION


vector<paraminit> initvals = {
    {"beta",1.96078,false},
    {"pdettheta",0.0503791,false},
    {"pdetcoef",2,true},
    {"imunity",1,true},
    {"shift",7.00153,false},
    {"mu",0.000575525,false},
    {"omega1",4.06449,false},
    {"omega2",0,true},
    {"rho",0.00249239,false},
    {"daktela",0.000449692,false},
    {"gammaa",0,true},
    {"gammas",0,true},
    {"gammad",0,true},
    {"varfactor",71.418,false},

};
#endif
#ifdef WORKING

static constexpr seirmodel::contrasttype ct = seirmodel::ectwls;  //// ectwls; // ;:;ectll;
static constexpr bool nopred = true;
static constexpr bool estimate = true;
static constexpr double timest = 10*60;
static constexpr bool sensanal = true;

static constexpr unsigned horizon=302;
static constexpr unsigned dusekpenalty=00;

#define DEADFROMIS
#define NOEXTRAMIMPORT
//#define PARKS
#define IMUNITY
#define DAKTELA
#define VARFACTOR
//#define POSITIVITY
#define EXPAPPROX
//#define LINAPPROX
//#define DONOTSUM
//#define QUARANTINERATIO
//#define VACATIONDUMMY
//#define TEMPER
#define SHIFT
//#define FACTORPAQ
#define MU
//#define FAMILYCONTACTS
//#define BETAADD
//#define ARATE
#define ASYMP
//#define CONSTIOTA
//#define PIECEWISEIOTA
//#define SETA
//#define PIECEWISEETA
//#define GAMMAS

//#define PIECEWISETHETA
//#define THREEPIECETHETA
//#define SCURVE

//#define K
//#define KACTIVE
//#define K2
#define WEIGHTVACCATION

vector<paraminit> initvals = {
    {"beta",0.669447,false},
    {"imunity",1,true},
    {"shift",0,true},
    {"mu",0.000322727,false},
    {"omega1",1.85324,false},
    {"omega2",0,true},
    {"rho",0.0182556,false},
    {"eta",0.1,true},
    {"daktela",1.89808e-37,true},
    {"varfactor",63.218,false},
};



#endif


#ifdef TR2V1

static constexpr seirmodel::contrasttype ct = seirmodel::ectll; // ectwls; // ;:;
static constexpr bool nopred = true;
static constexpr bool estimate = false;
static constexpr double timest = 4*60*60;

static constexpr unsigned horizon=246;
static constexpr unsigned dusekpenalty=00;

//#define VARFACTOR
//#define POSITIVITY
//#define EXPAPPROX
#define LINAPPROX
#define DONOTSUM
//#define QUARANTINERATIO
//#define VACATIONDUMMY
//#define TEMPER
//#define SHIFT
//#define FACTORPAQ
//#define MU
#define FAMILYCONTACTS
//#define BETAADD
#define ARATE
#define ASYMP
//#define CONSTIOTA
//#define PIECEWISEIOTA
//#define SETA
//#define PIECEWISEETA
//#define GAMMAS

#define PIECEWISETHETA
//#define SCURVE

//#define K
//#define KACTIVE
//#define K2


vector<unsigned> thetafrns({20,60,128,160,200,240});

vector<paraminit> initvals = {
{"marchimportrate",9.96795,true},// (0.801187)
{"beta",0.618509,false},// (0.0352657)
{"asymprate",0.549047,false},// (0.0256539)
{"omega",0.550946,false},// (0.0613134)
{"rho0",0.0194497,true},// (0.0147655)
{"rho1",0.00402584,true},// (0.00305399)
{"rho2",0.0204434,true},// (0.00610917)
{"rho3",0.0260418,true},// (0.00449858)
{"rho4",0.0248783,true},// (0.00366737)
{"rho5",0.00966628,true},// (0.000578787)
{"rho6",0.0221601,true},// (0.00148755)
{"eta",0.624971,true},// (0.00699778)
};



#endif
#ifdef TRV1

static constexpr seirmodel::contrasttype ct = seirmodel::ectwls; // ectwls; // ;:;
static constexpr bool nopred = true;
static constexpr bool estimate = false;
static constexpr double timest = 30*60;

static constexpr unsigned horizon=252;
static constexpr unsigned dusekpenalty=00;


//#define VARFACTOR
//#define POSITIVITY
#define EXPAPPROX // V2
// V2 #define LINAPPROX
#define DONOTSUM
//#define QUARANTINERATIO
//#define VACATIONDUMMY
//#define TEMPER
#define SHIFT
//#define FACTORPAQ
#define MU
//#define FAMILYCONTACTS
//#define BETAADD
//#define ARATE
#define ASYMP
//#define CONSTIOTA
//#define PIECEWISEIOTA
//#define SETA
//#define PIECEWISEETA
//#define GAMMAS

//#define PIECEWISETHETA
//#define SCURVE

//#define K
//#define KACTIVE
//#define K2




vector<paraminit> initvals = {
    {"marchimportrate",1.82882,false},// (0.204276)
    {"beta",1.19497,false},// (0.302035)
    {"shift",5.14679,false},// (0.314026)
    {"mu",0.00198165,false},// (0.000346351)
    {"omega1",2.30022,false},// (1.00413)
    {"omega2",0,true},// (1.01918)
    {"rho",0.0347199,false},// (0.000134977)
    {"eta",0.644159,false},// (0.00452093)
};

#endif

enum eparams {
#if !defined(NOEXTRAMIMPORT)
    emarchimportrate,
#endif
               ebeta,
#ifdef BETAADD
    ebetaadd,
#endif
#ifdef PARKS
    eparks,
#endif
#ifdef PDETTHETA
    epdettheta,
#endif
#ifdef PDETETA
    epdeteta,
#endif
#ifdef PDETCOEF
    epdetcoef,
#endif
#ifdef IMUNITY
    eimunity,
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
#ifdef LINAPPROX
    eomega,
#else
    eomega1,  eomega2,
#endif
#endif
#endif // PIECEWISEIOTA

#ifdef PIECEWISETHETA
               erho0,erho1, erho2, erho3, erho4, erho5,
               erho6, elastrho = erho6,
#else
#ifdef THREEPIECETHETA
    erho0,erho1, erho2, erho3, elastrho = erho3,
#else
#ifdef SCURVE
        erhoinit, erhosize, erhomid, erhok,
#else
        erho,
#endif
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
#if !defined(PDET)
               eeta,
#endif
#endif
#endif
#ifdef DAKTELA
    edaktela,
#endif
#ifdef K
    ek,
#endif

#ifdef K2
    ek2,
#endif

#ifdef GAMMAS
               egammaa,egammas,egammad,
#endif
#ifdef VARFACTOR
      evarfactor,
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



#if defined(PIECEWISETHETA) || defined(THREEPIECETHETA)
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
    enum eexp { eI, econtacts, egammared, ecalls, epositivity, egparks, epdet, edayadjust,
                enumexp};
private:
    static dmatrix getH()
    {
        dmatrix ret(czseir::enumobs,eseirnumstates);
        ret.setZero();
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

    dmatrix P(unsigned  t, const vector<double>& params, const G& g) const
    {
        double asymprate;
        expparams p;
#ifdef ARATE
        asymprate = params[easymprate];
#else
        asymprate = p.asymptomatic_rate;
#endif
        double sigma = p.sigma; //ksigma];
        dmatrix ret(eseirnumstates,eseirnumstates);
        ret.setZero();
        ret(eseire,eseirin) = sigma * asymprate;
        ret(eseire,eseiria) = sigma * (1-asymprate);
        double thetas = min(1.0,getthetas(t,params,g) * g.z[t][czseir::edayadjust]);
        double thetaa = getthetaa(t,params,g) * g.z[t][czseir::edayadjust];
        ret(eseire,eseirqe) = thetaa;
        ret(eseirqe,eseirinda) = sigma * asymprate;
        ret(eseirqe,eseiriada) = sigma * (1-asymprate);

        ret(eseirin,eseirjn) = p.delta_n;
        ret(eseirjn,eseirr) = p.gamma_In;
        ret(eseiria,eseiris) =
                p.symptoms_manifest_rate;

        ret(eseirjs,eseirr) = p.gamma_Is;
#ifdef MU
        ret(eseirjs,eseird) = params[emu];
#else
        ret(eseirjs,eseird) = p.mu;
#endif
        ret(eseiria,eseiriada) = thetaa;
        ret(eseirin,eseirinda) = thetaa;

        double mu;
#ifdef MU
        mu = params[emu];
#else
        mu = p.mu;
#endif

double muinis;
#ifdef DEADFROMIS
        ret(eseirisds,eseirdds) = 2*mu;
        ret(eseirisda,eseirdda) = 2*mu;
        muinis = 2*mu;
#else
        muinis = 0;
#endif
        if(thetas + muinis
                + p.delta_s <= 1)
        {
            ret(eseiris,eseirisds) = thetas;
            ret(eseiris,eseirdds) = muinis;
            ret(eseiris,eseirjs) = p.delta_s;
        }
        else if(thetas > 1) // nemelo by nastata
        {
            ret(eseiris,eseirisds) = 1,
            ret(eseiris,eseirjs) = 0;
        }
        else if(thetas + muinis>1)
        {
            ret(eseiris,eseirisds) = thetas;
            ret(eseiris,eseirdds) = 1-thetas;
            ret(eseiris,eseirjs) = 0;
        }
        else
        {
            ret(eseiris,eseirisds) = thetas,
            ret(eseiris,eseirdds) = muinis;
            ret(eseiris,eseirjs) = 1-thetas-muinis;;
        }

        ret(eseirinda,eseirjnda) = p.delta_n;
        ret(eseirjnda,eseirrda) = p.gamma_In;
        ret(eseiriada,eseirisda) = p.symptoms_manifest_rate;
        ret(eseirisda,eseirjsda) = p.delta_s;
        ret(eseirisds,eseirjsds) = p.delta_s;
        ret(eseirjsda,eseirrda) = p.gamma_Is;
        ret(eseirjsds,eseirrds) = p.gamma_Is;

//#if(!defined(DEADFROMIS))
        ret(eseirjsda,eseirdda) = mu;
        ret(eseirjsds,eseirdds) = mu;
        ret(eseirjs,eseirdds) = mu;
//#endif

        for(unsigned i=0; i<eseirnumstates; i++)
        {
            double s=0;
            for(unsigned j=0; j<eseirnumstates; j++)
                if(i != j)
                    s+=ret(i,j);
            if(s>1)
            {
                if(s > 1)
                {
                    cerr << m2csv(ret) << endl;
                    throw orpp::exception("s > 1!");
                }
                ret(i,i) = 0;
            }
            else
                ret(i,i) = 1.0-s;
        }


        return ret;
    }
    vector<double> I(unsigned  t, const vector<double>& params, const G& g) const
    {
        vector<double> ret(eseirnumstates,0);
#ifdef NOEXTRAMIMPORT
        double coef = 1.0 / (1.0 - 0.179);
#else
        double coef = t < 30 ? params[emarchimportrate] : 1.0 / (1.0 - 0.179);
#endif
        ret[eseire]= coef * g.z[t][eI];
        return ret;
    }


private:
    virtual vector<double> getiota(unsigned t, const vector<double>& p, const G& g) const = 0;
    virtual double getthetas(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual double getthetaa(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual string statelabel(unsigned i) const { return seirstatelabel(i); }
public:
    double numagpositive(const dvector& v) const
    {
        double s=0;

        for(unsigned i=0; i< eseirnumstates; i++)
        {
            switch(i)
            {
            case eseiriada:
            case eseirinda:
            case eseirisda:
            case eseiria:
            case eseiris:
            case eseirin:
            case eseirisds:
                s += v[i];
                break;

            case eseirjs:
            case eseirjn:
            case eseirjsda:
            case eseirjnda:
            case eseirjsds:
                s += v[i] / 2.0;
                break;

            }

        }
        return s;
    }

    double numantibodiesv(const dvector& v) const
    {
        double s=0;

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
                    s += v[i];
            }
        }
        return s;
    }

    double numagpositive(const evalresult& r, unsigned t) const
    {
        return numagpositive(r.est[t].x());
    }


    double numantibodies(const evalresult& r, unsigned t) const
    {
        return numantibodiesv(r.est[t].x());
    }

    double numactive(const G& g, unsigned t) const
    {
        double s=0;
        dvector est = g.est[t];
        for(unsigned i=0; i<= eseirin ; i++)
            s += est[i];
        return s;
    }

#ifdef VARFACTOR
    virtual double varfactor(const vector<double>& params ) const
    {
        return params[evarfactor];
    }
#endif

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


    dmatrix J(unsigned t, const vector<double>& params, const G& g) const
    {
        vector<double> iot = getiota(t,params,g);

        dmatrix ret(k(),k());
        ret.setZero();
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
    virtual dmatrix Gamma(unsigned /* t */, const vector<double>& params, const G& ) const
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

        dmatrix ret(m(),k());
        ret.setZero();
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
#ifdef GENERALTEST
    virtual void adjust(unsigned t, dvector& x) const
    {
        double attendance = 1;
        double sensitivity = 0.5;
        if(t==265  || t==264+14 )
        {
            double iap = sensitivity * attendance * x[eseiria];
            double inp = sensitivity *attendance * x[eseirin];
            double isp = sensitivity *attendance * x[eseiris];
            x[eseiria] -= iap;
            x[eseiriada] += iap;
            x[eseirin] -= inp;
            x[eseirinda] += inp;
            x[eseiris] -= isp;
            x[eseirisds] += isp;
            cout << iap + inp + isp << " removed from system" << endl;
        }
    }
#endif

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

#ifdef PARKS
       ba -=  p[eparks] * p[egparks];
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
#else // PIECEWISWIOTA
  #ifdef FACTORPAQ
double preiota = // exp(-p[eomega0] - p[eomegag]*g.z[t1][egammared] - p[eomegas]*g.z[t1][esred]);
                zofunction(1.0-p[eomegas]*(w* g.z[t1][esred] + (1-w)* g.z[t2][esred]))
                * zofunction(p[eomega0] -p[eomegag]*(w* g.z[t1][egammared] + (1-w)* g.z[t2][egammared]));
  #else // FACTORPAQ
        double c = w*g.z[t1][econtacts]+(1-w)*g.z[t2][econtacts];
     #ifdef FAMILYCONTACTS
        c -= p[efcontacts];
     #endif
     #ifdef VACATIONDUMMY
        if(t1>128 && t1 < 189)
            c *= p[evdummy];
     #endif

        double gammared = w*g.z[t1][egammared]+(1-w)*g.z[t2][egammared];
     #ifdef LINAPPROX
     #else
       #ifdef EXPAPPROX
        double preiota =c * exp(-p[eomega1]*gammared - p[eomega2]*gammared*gammared);
       #else
        double preiota =c * pow(1-p[eomega1]*gammared, p[eomega2]);
       #endif
        if(preiota  > 1.5)
            throw("preiota");
     #endif // LINAPROX
    #endif // FACTORPAQ
  #endif // PIECEWISEIOTA
  #ifdef CONSTIOTA
        preiota += p[eiotaconst];
  #endif
        double d;
  #ifdef TEMPER
        const unsigned t0 = -24;
        d = (1 - p[ephi]) + p[ephi]* cos((at-t0) / 365.0 * 2 * 3.141592 ) / cos((-t0) / 365.0 * 2 * 3.141592 ) ;
  #else
        d = 1;
#endif
        vector<double> iota(eseirnumstates,0);

#ifdef LINAPPROX
        double iot = max(0.0,c * ((p[ebeta]+ba) -p[eomega]*gammared )*d);
#else
        assert(preiota >= 0);
        double iot = max(p[ebeta]+ba,0.0)*d * preiota;

#ifdef IMUNITY
        iot *= (1-p[eimunity] * max(numantibodiesv(g.est[at]) / 10693000.00,0.0));
/*clog << at << "->" << p[eimunity] << "*"
     << numantibodiesv(g.est[at]) / 10693000.00 << "="
     << p[eimunity] * numantibodiesv(g.est[at]) / 10693000.00
     << " iot=" << iot << " c=" << c
      << " gammared=" << gammared << " e=" << exp(-p[eomega1]*gammared)
     << endl; */
#endif

#endif

        iota[eseiria] = iota[eseiris] = iot;
        iota[eseirin] = iot;
        assert(iota[eseiria] >= 0);
        assert(iota[eseirin] >= 0);
        assert(iota[eseiris] >= 0);
        return iota;
    }

    double weekave(unsigned t, const history& g) const
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

    double weekaved(unsigned t, const history& g) const
    {
        unsigned former = t>7 ? t-7 : 0;
        double numrep;
        if(t == 0)
            numrep=0;
        else
        {
            numrep = (g.y[t][eD]-g.y[former][eD])/ (t-former);
        }
        return numrep;
    }

    virtual double getthetas(unsigned t, const vector<double>& params, const G& g) const
    {        unsigned t1, t2;
             double w;
     #ifdef SHIFT
             double tt = t-params[eshift];
             t1 = max(0.0,floor(tt));
             t2 = min(t1+1,t);
             w = t2-tt;
     #else
             t1=t2=t;
             w=1;
     #endif


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
#ifdef PDET
        double eta = 0;
#else
        double eta = params[eeta];
#endif
#endif
#ifdef DONOTSUM
        rh = eta;
#else
        rh = getthetaa(t,params,g) + eta;
#endif
#endif
#ifdef POSITIVITY
      double positivity = g.z[t][epositivity];
      rh *= 1-params[eetaposcoef] * positivity;
#endif
       expparams p;
       double varrhos = p.delta_s + params[emu];
       double pdet = g.z[t][epdet];
       assert(pdet <=1.0);
#ifdef PDETETA
       rh += params[epdeteta] * pdet;
#endif
#ifdef PDET

#ifdef PDETCOEF
                 pdet*=params[epdetcoef];
#endif
       assert(pdet <=1.0);
       assert(pdet >= 0.0);


#ifdef PDETONLYETA
       rh = pdet * varrhos / (1-pdet);
#else
        double varrhos = p.delta_s + params[emu];
        double varrhoa = p.delta_n;
        double theta = getthetaa(t,params,g);
        double alpha = p.asymptomatic_rate;
        double varsigma = p.symptoms_manifest_rate;
        double p1 = alpha * theta / (theta + varrhoa);
        double p2 = (1-alpha) * theta / (theta+ varsigma);
        double q =  (1-´alpha) * varsigma / (theta + varsigma);
        double x = q / (pdet
                        - p1 - p2);
        rh = max(0.0, varrhos / x);
#endif
#endif
        assert(rh >= 0);
        return zofunction(rh,0.000001);
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
#if defined(PIECEWISETHETA) || defined(THREEPIECETHETA)

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
      double positivity = g.z[t][epositivity];
      c *= 1-params[ethetaposcoef] * positivity;
#endif

#ifdef DAKTELA
      double wa = weekave(t,g);
      c += params[edaktela] * (wa > 0 ? g.z[t][ecalls] / wa : 0);
#endif
#ifdef PDETTHETA
      c += params[epdettheta] * g.z[t][epdet];
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
            si.y.push_back(dv(obs));
        }
        else
            rfin = true;
        vector<double> exp;
        for(unsigned j=0; j<czseir::enumexp; j++)
            exp.push_back(c.getdouble(i+1,1+firstexp+j));
        si.z.push_back(dv(exp));
    }
}



void seir()
{
    csv<','> c("../input/cz.csv");

    obsseir s;

    seirmodel::data si;

    fillsi(si,c,0);

    for(unsigned i=0; i<si.y.size(); i++)
    {

        double a = s.weekave(i,si);
        double d = s.weekaved(i,si);
        double w, wd;
        if(a < 20)
            w =1.0 / 20.0;
        else
            w=1 / a;
        if(d < 10)
            wd = 1 / 10.0;
        else
            wd = 1 / d;
#ifdef WEIGHTVACCATION
        if(i < 180)
        {
            w /=20;
            wd /= 20;
        }
#endif
        si.w.push_back({w,w,wd});
    }

    assert(horizon <= si.y.size());
    si.y.resize(horizon);
    si.z.resize(min(100UL+horizon,si.z.size()));

    auto siest = si;





    vector<paraminfo> params
            = {
#if !defined(NOEXTRAMIMPORT)
                 paraminfo("marchimportrate",3,0,10),
#endif
                 paraminfo("beta",2.08,0,8),
#ifdef BETAADD
              paraminfo("betaadd",0,0,3),
#endif
#ifdef PARKS
              paraminfo("parks",0,-1,1),
#endif
#ifdef PDETTHETA
        paraminfo("pdettheta",0,0,1),
#endif
#ifdef PDETETA
        paraminfo("pdeteta",0,0,1),
#endif
#ifdef PDETCOEF
        paraminfo("pdetcoef",1,0,2.5),
#endif

#ifdef IMUNITY
             paraminfo("imunity",0,0,10),
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
            paraminfo("shift",5,0,15),
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
    #ifdef LINAPPROX
        paraminfo("omega", 0.05, 0, 1),
    #else
                paraminfo("omega1", 2, 0, 10),
                paraminfo("omega2", 0, 0, 10),
    #endif
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
#ifdef THREEPIECETHETA
        paraminfo("rho0", 0, 0, 0.7),
        paraminfo("rho1", 0.1, 0, 0.7),
        paraminfo("rho2", 0.1, 0, 0.7),
        paraminfo("rho3", 0.2, 0, 0.7),
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
#endif

#ifdef POSITIVITY
        paraminfo("thetaposcoef", 0.5, 0, 1),
        paraminfo("etaposcoef", 0.5, 0, 1),
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
    #ifdef PDET
    #else
        paraminfo("eta", 0.3, 0, 0.9),
    #endif
#endif
#endif
#ifdef DAKTELA
        paraminfo("daktela", 0.025, 0, 0.1),
#endif

#ifdef K
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
                paraminfo("gammaa",0, 0 , 0.0001 ),
                paraminfo("gammas",0, 0 , 0.0001 ),
                paraminfo("gammad",0,0, 100),
#endif
#ifdef VARFACTOR
        paraminfo("varfactor",1,1,150),
#endif

              };

assert(params.size()==enumparams);

    vector<bool> filter(params.size(),false);
//    filter[eomega0]=filter[eomega1]=filter[efcontacts]=true;
//    for(unsigned i=erho0; i<=elastrho; i++)
//        filter[i] = true;
//    for(unsigned i=erho0; i<=elastrho; i++)
//        filter[i] = true;


    vector<double> rp(enumparams);

    for(unsigned i=0; i<params.size(); i++)
    {
        rp[i] = params[i].initial;
        for(unsigned j=0; j<initvals.size(); j++)
            if(params[i].name == initvals[j].n)
            {
                rp[i]=params[i].initial = initvals[j].v;
                filter[i] = initvals[j].filtered;
                clog << params[i].name << " set filter to " << filter[i] << endl;
            }
    }


    if(1)
    {
        if(estimate)
        {
            uncertain res;
            clog << "ll= " << s.estimate(params,siest,res,ct,filter,timest) << endl;
            rp = vd(res.x());
            vector<double> g = s.grad(rp,si,ct);
            for(unsigned i=0; i<params.size(); i++)
            {
                clog << "{\"" << params[i].name << "\","
                    << rp[i] << "," << (filter[i] ? "true" : "false") << "},";
                clog << "// (" << sqrt(res.var()(i,i)) << ") " << g[i] << endl;
            }
            clog << endl << "var of params" << endl;
            clog << res.var() << endl;
        }

        clog << dv(rp).transpose() << endl;

        czseir::evalresult r = s.eval<true>(rp,si, si.z.size()-si.y.size(),ct);
        ofstream o("output.csv");
        r.output(o,si);
        clog <<  "ll=" << r.contrast  << endl
             <<  "ll/n=" << r.contrast / r.contrasts.size() << endl
             << "additional contrast = " << s.contrastaddition(r) << endl
             << "Duseks survey = " << s.numantibodies(r,70) << endl
             << "Teachars survey = " << s.numagpositive(r,280) << endl
             << "People with antibodies = " << s.numantibodies(r,289) << endl
             << endl << endl;

        if(0) // difference
        {
            auto rrr = rp;
            for(unsigned i=0; i<20; i++)
            {
               czseir::evalresult rr = s.eval<true>(rrr,si, si.z.size()-si.y.size(),ct);
               cout << rrr[erho] << "->" << rr.contrast << endl;
               rrr[erho] += 0.000001;
            }
            throw;
//            rrr[epdeteta] += 0.000001;
            czseir::evalresult rr = s.eval<true>(rrr,si, si.z.size()-si.y.size(),ct);
            for(unsigned i=0; i<rr.contrasts.size(); i++)
            {
                clog << i << ","
                     << r.contrasts[i] << "," << rr.contrasts[i]
                     << "," << r.contrasts[i]-rr.contrasts[i] << endl;
            }
            throw;
        }
        if(sensanal)
        {
            clog << "Sensitivity analysis" << endl;
            auto r = s.sensitivity(rp,si,ct);
            for(unsigned i=0; i<r.first.size(); i++)
            {
                double m = r.second[i][r.second[i].size()/2];
                clog << params[i].name;
                unsigned h = r.first[i].size() / 2;
                double s1=0;
                double s2=0;
                for(unsigned j=0; j<r.first[i].size(); j++)
                {
                    if(j<h)
                        s1+=r.second[i][j]-m;
                    else if(j>h)
                        s2+=r.second[i][j]-m;
                    clog << "," << r.second[i][j]-m;
                }
                clog << ",," << r.first[i][0] << ","
                     << r.first[i][r.first[i].size()-1] << ",";
                if(s1<0 && s2<0)
                    clog << "OK";
                else if(s1>0)
                    clog << "<-";
                else
                    clog << "->";
                clog << endl;
            }
        }
        if(0) // confidence of reported
        {
            cout << "stddevs R" << endl;
            for(unsigned i=0; i<r.pred.size(); i++)
            {
                dmatrix v = r.pred[i].var();
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
            dmatrix T(eseirin+1,eseirin+1);
            T.setZero();
            for(unsigned i=0; i<7; i++)
                T = T + r.Ts[si.y.size()-1-i].block(0,0,eseirin+1,eseirin+1);
            T = (1.0 / 7.0) * T;
            double gamma = T(0,1);
            double c=si.z[si.z.size()-1][czseir::econtacts];
            cout << "C,rho" << endl;
            for(unsigned i=30; i<100; i++)
            {
                double factor = (i/100.0) / c;
                T(0,1)=T(0,2)=T(0,3)=gamma*factor;
                cout << i/100.0 << "," << radius(T) << endl;
            }
            throw;
        }



        clog <<
#ifdef ASYMP
                "RA,predRA,predRAtrend, resA, resAtrend,"
                "RS,predPS, predRStrend, resS, resStrend,"
#else
                "R,predR, predRtrend, res, restrend,,,,,,"
#endif
                "D,predD, predDtrend, resD, resDtrend, "
                "gamma, eta, theta, rho, rhotrend,"
                "R, Rcp ,contrast" << endl;
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
            double waminus7;
            double wadminus7;
            dmatrix trendm;
            if(i>15 && !nopred)
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
                trendm = rx.Ts[i].block(0,0,eseirin+1,eseirin+1);
                waminus7 = max(s.weekave(i-ds,six),1.0);
                wadminus7 = max(s.weekaved(i-ds,six),1.0);
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
                trendm = r.Ts[i].block(0,0,eseirin+1,eseirin+1);
                waminus7 = 1;
                wadminus7 = 1;
            }
#ifdef ASYMP
            double predPA = r.pred[i].x()[eseirnumstates+czseir::eRA];
            double predPS = r.pred[i].x()[eseirnumstates+czseir::eRS];
#else
            double predP = r.pred[i].x()[eseirnumstates+czseir::eR];
#endif
            double predD = r.pred[i].x()[eseirnumstates+czseir::eD];

            double rng = radius(r.Ts[i].block(0,0,eseirin+1,eseirin+1));

            double rn = s.R(i,r);
            double rcp = s.Rcp(i,r,4);

            clog
#ifdef ASYMP
                    << PA << "," << predPA << "," << predPAtrend << ","
                    << (PA - predPA)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eRA,eseirnumstates+czseir::eRA)) << ","
                    << (PA - predPAtrend ) / waminus7 << ","
                    << PS << "," << predPS << "," << predPStrend << ","
                    << (PS - predPS)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eRS,eseirnumstates+czseir::eRS)) << ","
                    << (PS - predPStrend) / waminus7 << ","
#else
                 << P << "," << predP << "," << predPtrend << ","
                 << (P-predP)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eR,eseirnumstates+czseir::eR)) << ","
                 << P - predPtrend  << ",,,,,,"
#endif
                 << D << "," << predD << "," << predDtrend << ","
                         << (D - predD) / sqrt(r.pred[i].var()(eseirnumstates+czseir::eD,eseirnumstates+czseir::eD)) << ","
                         << (D - predDtrend) / wadminus7 << ","
                 << r.Ts[i](eseire,eseiris) << ","
                 << r.Ts[i](eseirisds,eseiris) << ","
                 << r.Ts[i](eseirqe,eseire) << ","
                 << rng << ","
                << radius(trendm) << ","
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
    throw;

/*    auto old_clogbuf = std::clog.rdbuf();

    for(int ar=2; ar < 8; ar++)
    {
        cout << "starting " << ar << endl;
        ostringstream os;
        os << "log" << ar << ".txt";
        ofstream o(os.str());


        std::clog.rdbuf(o.rdbuf());

        try
        {
            seir(ar / 10.0);
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
    }
    clog.rdbuf(old_clogbuf);
    return 0;*/
}

