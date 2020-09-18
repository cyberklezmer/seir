
#include <vector>
#include "emodel.hpp"
#include "orpp.hpp"
#include "seir.hpp"
#include <cmath>

using namespace std;
using namespace orpp;

#define FAMILYCONTACTS
//#define SETA
//#define APAR
//#define ASYMP
#define GAMMAS
//#define PIECEWISE
//#define BETAADD
#define SCURVE
#define K
#define K2
//#define TESTSPERCAND

enum eparams {
               ebeta,
#ifdef BETAADD
    ebetaadd,
#endif
#ifdef FAMILYCONTACTS
    efcontacts,
#endif
    eomega, // ebetaadd,
#ifdef GAMMAS
               egammad, egammar,
#endif
#ifdef K
    ek,
#endif
#ifdef K2
    ek2,
#endif
#ifdef APAR
    easymprate,
#endif
#ifdef PIECEWISE
               erho0,erho1, erho2, erho3, erho4, erho5,
               erho6, elastrho = erho6,
#endif
#ifdef ETAPIECE
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

#ifdef SCURVE
//               erhomin,
    erhosize, erhomid, erhok,
#endif
#ifdef TESTSPERCAND
    erhocoef,
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

inline string seirstatelabel(unsigned i)
{
    static const string l[] = {
                   "E",
                   "Ia", "Is", "In", "Js", "Jn", "R", "D",
                   "QE",

                   "Iada", "Inda",
                   "Isda", "Jsda","Jnda", "Rda", "Dda"

                   "Isds", "Jsds","Rds", "Dds"
    };
    return l[i];
}


class czseir : public seirmodel
{
public:
    enum eobs {
#ifdef ASYMP
        eRA, eRS,
#else
        eT,
#endif
        eD,  enumobs};
    enum eexp { eI, econtacts, emasks, ehygiene, edayadjust, etests, esmT,
                enumexp};
private:
    static matrix getH()
    {
        matrix ret(czseir::enumobs,eseirnumstates,0);
#ifdef ASYMP
        for(unsigned i=eseiriada; i<=eseirdda; i++)
            ret(eRA,i) = 1;
        for(unsigned i=eseiriada; i<=eseirdds; i++)
            ret(eRS,i) = 1;
        ret(eRS,eseird)=1;
#else
        for(unsigned i=eseiriada; i<eseirnumstates; i++)
            ret(eT,i) = 1;
        ret(eT,eseird)=1;
#endif
        ret(eD,eseird)=1;
        ret(eD,eseirdda)=1;
        ret(eD,eseirdds)=1;
        return ret;
    }
    static matrix getU()
    {
        matrix ret(eseirnumstates,enumexp,0);
        expparams p;

        ret(eseire,eI)=1.0 / (1.0 - p.asymptomatic_rate);
        return ret;
    }
public:
    czseir() : seirmodel(getH(),getU()) {}

    matrix P(unsigned  t, const vector<double>& params, const G& g) const
    {
        double asymprate;
        expparams p;
#ifdef ARATE
        asymprate = pararams[easymprate];
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
        ret[eseirjs][eseird] = p.mu;

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
        ret[eseirjsda][eseirdda] = p.mu;
        ret[eseirjsds][eseirdds] = p.mu;

        for(unsigned i=0; i<eseirnumstates; i++)
        {
            double s=0;
            for(unsigned j=0; j<eseirnumstates; j++)
                if(i != j)
                    s+=ret[i][j];
            if(s>1)
            {
                if(s > 1.000000001)
                    throw "s > 1!";
                ret[i][i] = 0;
            }
            else
                ret[i][i] = 1.0-s;
        }
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
        return ll;
    }


    matrix I(unsigned t, const vector<double>& params, const G& g) const
    {
        vector<double> iot = getiota(t,params,g);
        matrix ret(k(),k(),0);
        for(unsigned i=0; i<eseirnumstates; i++)
            ret(0,i) = iot[i];
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
        double gammar;
#ifdef GAMMAS
        gammad = params[egammad];
        gammar = params[egammar];
#else
        gammad = 0.01;
        gammar = 0.0000001;
#endif

        matrix ret = matrix(m(),k(),0);
#ifdef ASYMP
        ret(eRS,eseiris) = gammar;
        ret(eRA,eseire+eseiria) = gammar;
#else
        ret(eT,eseiris) = gammar;
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

    vector<double> getiota(unsigned t, const vector<double>& p, const G& g) const
    {
        assert(p.size()==enumparams);
        double ba=0;
#ifdef BETAADD
        ba = p[ebetaadd];
#endif
        double c = g.z[t][econtacts];
#ifdef FAMILYCONTACTS
        c -= p[efcontacts];
#endif
        double preiota =c*(1-0.85*g.z[t][emasks])*(1-p[eomega] *g.z[t][ehygiene]);
        vector<double> iota(eseirnumstates,0);
        iota[eseiria] = iota[eseiris] = (p[ebeta]/*+p[ebetaadd]*/)*preiota;
        iota[eseirin] = (p[ebeta]+ba)*preiota;
        return iota;
    }

    virtual double getthetas(unsigned t, const vector<double>& params, const G& g) const
    {
        double rh;
        double numrep = g.z[t][esmT];
#ifdef ETAPIECE
        unsigned numetas = elasteta + 1 - eeta0;
        unsigned period = 35;
        unsigned i = min(t / period, numetas - 1);
        if(i == numetas - 1)
            rh = params[eeta0 + i];
        else
        {
            double w = (t % period) / static_cast<double>(period);
            rh = (1-w) * params[eeta0+i] + w * params[eeta0+i+1];
        }
#else
#ifdef SETA
        double eta = params[eetamin]
                + params[eetasize] / (1.0 + exp(-params[eetak] * (t -params[eetamid] )));
#else
        double eta = params[eeta];
#endif
        rh = getthetaa(t,params,g) + eta;
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
        double c = numrep < k ? rh : rh * k / numrep ;
if(c> 0.7)
    cout << "c > 0.7";
        return max(0.0,min(0.7,c));
    }


    virtual double getthetaa(unsigned t, const vector<double>& params, const G& g) const
    {
#ifdef TESTSPERCAND
        double numtraceable = g.est[t][eseire]+g.est[t][eseiria]
                                +g.est[t][eseiris]+g.est[t][eseirin];
        if(numtraceable == 0)
            return 0;
        else
            return min(params[erhocoef] * g.z[t][etests]/numtraceable,0.7);
#endif
        double rh;
        double numrep = g.z[t][esmT];
#ifdef PIECEWISE
        unsigned numrhos = elastrho + 1 - erho0;
        unsigned period = 35;
        unsigned i = min(t / period ,numrhos - 1);
        if(i == numrhos - 1)
            rh = params[erho0 + i];
        else
        {
            double w = (t % period) / static_cast<double>(period);
            rh = (1-w) * params[erho0+i] + w * params[erho0+i+1];
        }
#endif
#ifdef SCURVE
       rh = //params[erhomin]
               params[erhosize] / (1.0 + exp(-params[erhok] * (t -params[erhomid] )));
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
        double c = numrep < k ? rh : rh * k / numrep ;

if(c> 0.7)
   cout << "c > 0.7";
        return max(0.0,min(0.7,c));
    }
};

void fillsi(seirmodel::data& si, csv<','>c, unsigned lag = 0)
{
    bool rfin = false;
    for(unsigned j=0; j<czseir::enumobs; j++)
    {
        si.ylabels.push_back(c(0,1+j));
    }

    for(unsigned j=0; j<czseir::enumexp; j++)
        si.zlabels.push_back(c(0,czseir::enumobs+1+j));

    for(unsigned i=lag; i< c.r()-1; i++)
    {

        si.dates.push_back(c(i+1,0));
        unsigned k=1;
        if(c(i+1,1) != "" && !rfin)
        {
            vector<double> obs;
            for(unsigned j=0; j<czseir::enumobs; j++)
                obs.push_back(c.getunsigned(i+1,k++));
            si.y.push_back(obs);
        }
        else
        {
            rfin = true;
            k += czseir::enumobs;
        }
        vector<double> exp;
        for(unsigned j=0; j<czseir::enumexp; j++)
//if(j!=czseir::econtacts && j!=czseir::emasks)
            exp.push_back(c.getdouble(i+1,k++));
//else
//  exp.push_back(c.getdouble(i+10 < c.r() ? i+10 : c.r()-1,k++));
        si.z.push_back(exp);
    }
}



void seir()
{
#ifdef ASYMP
    csv<','> c("../input/cza.csv");
#else
    csv<','> c("../input/cz.csv");
#endif
    seirmodel::data si;

    fillsi(si,c,0);
//    si.y.resize(si.y.size()-70);
//    si.z.resize(si.z.size()-70);

    auto siest = si;

    vector<double> initvals = {
        2.66452, // beta
        0.146487, // fcontacts
        0.889302, // omega
        1.6263e-17, // gammad
        22.9477, // gammar
        105.594, // k
        131.547, // k2
        0.107862, // eta
        0.414669, // rhosize
        129.49, // rhomid
        0.0677453, // rhok
    };

    vector<paraminfo> params
            = {
                 paraminfo("beta",1.47992,0,3),
#ifdef BETAADD
              paraminfo("betaadd",0,0,3),
#endif
#ifdef FAMILYCONTACTS
                paraminfo("fcontacts",0.16,0,0.4),
#endif

                paraminfo("omega", 0.828601, 0.5, 1),
#ifdef GAMMAS
                paraminfo("gammad",0,0, 100),
                paraminfo("gammar",10, 0 , 100 ),
#endif
#ifdef K
                paraminfo("k", 50, 0, 260),
#endif
#ifdef K2
               paraminfo("k2", 50, 0, 325),
#endif
#ifdef APAR
        paraminfo("asymprate", 0.19, 0, 1),
#endif
#ifdef PIECEWISE
                paraminfo("rho0", 0, 0, 0.7),
                paraminfo("rho1", 0.05, 0, 0.7),
                paraminfo("rho2", 0.05, 0, 0.7),
                paraminfo("rho3", 0.06, 0, 0.7),
                paraminfo("rho4", 0.20, 0, 0.7),
                paraminfo("rho5", 0.50, 0, 0.7),
                paraminfo("rho6", 0.5, 0, 0.7),
#endif
#ifdef ETAPIECE
                paraminfo("eta0", 0, 0, 0.7),
                paraminfo("eta1", 0.05, 0, 0.7),
                paraminfo("eta2", 0.05, 0, 0.7),
                paraminfo("eta3", 0.06, 0, 0.7),
                paraminfo("eta4", 0.20, 0, 0.7),
                paraminfo("eta5", 0.50, 0, 0.7),
                paraminfo("eta6", 0.5, 0, 0.7),
#else
    #ifdef SETA
        paraminfo("etamin", 0, 0, 0.3),
        paraminfo("etasize", 0.0, -0.7, 0.7),
        paraminfo("etamid", 100, 0, 300),
        paraminfo("etak", 0.06, 0, 1),
    #else
               paraminfo("eta", 0, 0, 0.7),
#endif
#endif
#ifdef SCURVE
        paraminfo("rhosize", 0.3, 0, 0.7),
        paraminfo("rhomid", 100, 0, 300),
        paraminfo("rhok", 0.06, 0, 1),
#endif
#ifdef TESTSPERCAND
        paraminfo("rhocoef", 0.1, 0, 1)
#endif

              };
    vector<bool> filter(params.size(),false);
    obsseir s;

//    for(unsigned i=0; i<params.size(); i++)
//        params[i].initial = initvals[i];
    if(1)
    {
        vector<double> rp;
        if(1)
        {
            uncertain res;
            cout << "ll= " << s.estimate(params,siest,res,filter) << endl;
            for(unsigned i=0; i<params.size(); i++)
            {
                cout << params[i].name << " ";
                cout << res.x()[i] << "(" << sqrt(res.var()(i,i)) << ")" << endl;
                params[i].initial = res.x()[i];
                filter[i]= true;
            }
            rp = res.x();
        }
        else
            rp = initvals;

        csvout<double,','>(clog, rp);
        czseir::evalresult r = s.eval<true>(rp,si);
        ofstream o("output.csv");
        r.output(o,si);
        clog << "ll/n=" << r.contrast / r.contrasts.size() << endl
             << "additional contrast = " << s.contrastaddition(r) << endl
             << "Duseks survey = " << s.numantibodies(r,70) << endl
             << endl << endl;

        if(0) // zkroceni
        {
            matrix T(eseirin+1,eseirin+1,0);
            for(unsigned i=0; i<7; i++)
                T = T + r.Ts[si.y.size()-1-i].submatrix(0,0,eseirin+1,eseirin+1);
            T = (1.0 / 7.0) * T;

            clog << "T" << endl <<  T << endl;
            clog << "rho(T)=" << range(T) << endl << endl;


            double c=si.z[si.z.size()-1][czseir::econtacts];
            double c0=rp[efcontacts];
            clog << "c=" << c << endl;

            double delta = 0.01;
            clog << "theta,rho,iota98" << endl;

            matrix I(4,4,0);
            I(0,0)=I(1,1)=I(2,2)=I(3,3)=1;

            vector<double> m(4,0);
            m[0] = 30;

            double required = 0.98;
            double theta = 0;
            for(unsigned i=0; i<7; i++)
                theta += r.Ts[si.y.size()-1-i](eseirqe,eseire);
            theta /= 7.0;
            bool firsttime = true;


            for(; theta<0.30; theta+= delta)
            {

                double kappa=T(0,1);
                double origkappa = kappa;
                double factor = origkappa / (c-c0);

                if(range(T) > required)
                {
                    auto Tt =T;
                    for(unsigned i=0; i<100000; i++)
                    {
                        kappa -= 0.001;
                        Tt(0,1)=Tt(0,2)=Tt(0,3)=kappa;
                        double cc = kappa / factor + c0;
                        if(firsttime && cc - floor(cc*100)/100 < 0.005)
                            clog << "," << cc << "," << range(Tt) << endl;

                        if(range(Tt) < 0.98)
                            break;
                    }
                    clog << theta << "," << range(T) << ",";
                    clog << kappa / factor + c0;
                }

                auto bound = ((I - T).inverse())*m;
                double s=0;
                for(unsigned i=0; i<4; i++)
                    s+=bound(i,0);
                clog << "," << s << endl;

                T(0,0) -= delta;
                T(1,1) -= delta;
                T(2,2) -= delta;
                T(3,3) -= delta;
                firsttime = false;
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
                "D,predD, predDtrend,dres, drestrend, iota, theta, theteq,r,rtrend,"
                "diota, drho,contrast" << endl;
        for(unsigned i=0; i<si.y.size(); i++)
        {
            double D = si.y[i][czseir::eD];
#ifdef ASYMP
            double PA = si.y[i][czseir::eRA];
            double PS = si.y[i][czseir::eRS];
            double predPAtrend;
            double predPStrend;
#else
            double P = si.y[i][czseir::eT];
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
                predPtrend = rx.pred[i].x()[eseirnumstates+czseir::eT];
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
            double predP = r.pred[i].x()[eseirnumstates+czseir::eT];
#endif
            double predD = r.pred[i].x()[eseirnumstates+czseir::eD];
//            double predPtrend = r2.pred[i].x()[eseirnumstates+czseir::eT];
            double delta = 0.0001;

            double rng = range(r.Ts[i].submatrix(0,0,eseirin+1,eseirin+1));
            matrix Tiota = r.Ts[i].submatrix(0,0,eseirin+1,eseirin+1);
            Tiota(eseire,eseiria) += delta;
            Tiota(eseire,eseirin) += delta;
            Tiota(eseire,eseiris) += delta;
            matrix Trho = r.Ts[i].submatrix(0,0,eseirin+1,eseirin+1);
            Trho(eseire,eseire) -= delta;
            Trho(eseiria,eseiria) -= delta;
            Trho(eseirin,eseirin) -= delta;
            double diota = (range(Tiota)-rng)/delta;
            double drho = (range(Trho)-rng)/delta;

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
                 << (predP - P)/sqrt(r.pred[i].var()(eseirnumstates+czseir::eT,eseirnumstates+czseir::eT)) << ","
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
                << diota << "," << drho;
            if(i<r.contrasts.size())
                clog << "," << r.contrasts[i];
            clog << endl;
        }
        for(unsigned i=0; i<params.size(); i++)
        {
            cout << rp[i] << ", // ";
            cout << params[i].name << " ";
//            cout << "(" << sqrt(res.var()(i,i)) << ")"
//
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
//        seir();
//  closure();
//        graph();
        seir();

/*        if(1)
            graph();
        else
            epi();
*/        cout << time(0) - startt << " seconds" << endl;

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

