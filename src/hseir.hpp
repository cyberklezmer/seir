#ifndef HSEIR_HPP
#define HSEIR_HPP

#include "cohortseir.hpp"
#include "estimableseirfilter.hpp"

class hpartial: public partialseir
{
public:

    enum states { E,  Ia,  Ip,  Is, Iu, numactives, R=numactives,

                  Edelta, Iadelta, Ipdelta,  Isdelta,  Rdelta,
                  Hdelta, Rhdelta,
                  Ddelta, Dhdelta,

                  Isd, Rd,
                  Hd, Rhd,
                  Dd, Dhd,

                  numstates};

    enum params { alpha, sigma, varsigma, gammas, gammaa, iotas, mus,
                  theta, eta, pi,
                  gammah, muh,
                  musratio,
                  prebeta,
                  betafactor,
                  ufactor,
                  vfactor,
                  numparams
                };

    hpartial() : partialseir()
    {}
    virtual dmatrix P(const vector<double>& params) const
    {
        /// more intutive to construct transposition
        dmatrix Pt(k(),k());
        Pt.setZero();

        double mh= params[muh] ;
        double ms = params[mus] * params[musratio];

        double A = params[theta] *
            (1
             + params[alpha] * params[sigma] / (params[theta] + params[gammaa])
             + (1-params[alpha]) * params[sigma] / (params[varsigma] +params[theta]))
          / (params[theta]+params[sigma]);
        double B =
              (1-params[alpha]) * params[sigma] * params[varsigma] * params[eta]
                      / (params[gammas]+ms+params[gammas]+params[eta])
                      / (params[theta]+params[sigma])
                      / (params[varsigma]+params[theta]);
        double pu = max(min(1.0, (params[pi] - A ) / B),0.0);


        Pt(E,Ia) = params[sigma] * params[alpha];
        Pt(E,Ip) = params[sigma] * (1-params[alpha]);
        Pt(Ip,Is) = params[varsigma] * pu;

        Pt(Ip,Hd) = params[iotas];

        Pt(Ip,Iu) = params[varsigma] *(1-pu);

        Pt(Ia,R) = params[gammaa];
        Pt(Is,R) = params[gammas];
        Pt(Iu,R) = params[gammas];

        Pt(E,Edelta) = params[theta];
        Pt(Ia,Iadelta) = params[theta];
        Pt(Ip,Ipdelta) = params[theta];
        Pt(Is,Isd) = params[eta];
        Pt(Is,Hd) = params[iotas];
        Pt(Is,Dd) = ms;
        Pt(Iu,Hd) = params[ufactor] * params[iotas];
        Pt(Iu,Dd) = params[vfactor] * ms;

        Pt(Edelta,Iadelta) = params[sigma] * params[alpha];
        Pt(Edelta,Ipdelta) = params[sigma] * (1-params[alpha]);
        Pt(Ipdelta,Isdelta) = params[varsigma];
        Pt(Iadelta,Rdelta) = params[gammaa];
        Pt(Isdelta,Rdelta) = params[gammas];
        Pt(Isdelta,Ddelta) = ms;
        Pt(Isdelta,Hdelta) = params[iotas];
        Pt(Hdelta,Rhdelta) = params[gammah];
        Pt(Hdelta,Dhdelta) = mh;

        Pt(Isd,Rd) = params[gammas];
        Pt(Isd,Dd) = ms;
        Pt(Isd,Hd) = params[iotas];
        Pt(Hd,Rhd) = params[gammah];
        Pt(Hd,Dhd) = mh;

        for(unsigned i=0; i<k(); i++)
        {
            double s = 0;
            for(unsigned j=0; j<k(); j++)
                s += Pt(i,j);
            Pt(i,i)=max(0.0,1-s);
        }

        if(0)
        {
            for(unsigned i=0; i< k(); i++)
                for(unsigned j=0; j< k(); j++)
                    if(Pt(i,j) < 0 || Pt(i,j) > 1)
                    {
                        cerr << m2csv(Pt) << endl;
                        throw "not a probability matrix";
                    }
        }
        return Pt.transpose();
    }


    virtual dmatrix hatB(const vector<double>& params) const
    {
        dmatrix Bt(k(),k());
        Bt.setZero();
        double b = params[prebeta]*params[betafactor];
        Bt(Is,E) = b;
        Bt(Ip,E) = b;
        Bt(Ia,E) = b / 4;
        Bt(Iu,E) = b;

        return Bt.transpose();
    }

    std::string statelabel(unsigned i) const
    {
        static std::vector<std::string> r =
        { "E", "Ia", "Ip" ,"Is","Iu","R",
          "Edelta","Iadelta","Ipdelta","Isdelta","Rdelta",
          "Hdelta","Rhdelta",
          "Ddelta","Dhdelta",
          "Isd","Rd",
          "Hd","Rhd",
          "Dd","Dhd" };
        assert(i<r.size());
        return r[i];
    }

    virtual unsigned k() const
    {
        return numstates;
    }

    virtual double susceptible( double population, const dvector&  compartments ) const
    {
        double s=0;
        for(unsigned i=0; i<compartments.size(); i++)
            s += compartments[i];
        assert(s < population);
        return 1.0 - s / population;
    }

    virtual dmatrix activesubmatrix( const dmatrix& T) const
    {
        return T.block(0,0,numactives,numactives);
    }

};

class hcohortseir: public cohortseir<hpartial>
{
public:
    enum excolumns {I0,	I20,	I65,	I80,
                     DAYADJUST,PDET,REDUCTIONMEAN,BRIGITATTACK,FEAR,BETAFACTOR,BFACTOR,
                     F0, F20, F65, F80, FALL, S0, S20, S65, S80, SALL,
                     numexcolumns};

    enum computationparams {
               paqshift,
               varbfactor,
               ce,
                         firstdisp=ce,
               cd,
               cu,
               ch,
                        lastdisp=ch,
               rcoef,
               hcoef,
               gcoef,
               dcoef,
                        lastvar=dcoef,
               omega,
               omega2,
               pi,
               thetacoef,
               eta0,
               theta0,
               hybrigitefect,
               hobrigitefect,
               brigiteeffect,
               cycleeffect,
               numcomputationparams
    };

    enum commonparams
    {
        sigma = numcomputationparams,
        varsigma,
        ufactor,
        vfactor,
        gammas,
        gammaa,
        mus,
        numcommonandcompparams
    };

    static std::vector<hpartial::params> commonpars()
    {
           return {hpartial::sigma,
               hpartial::varsigma,
               hpartial::ufactor,
               hpartial::vfactor,
               hpartial::gammas,
               hpartial::gammaa,
               hpartial::mus
              };
    }


    enum exclusiveparams
    {
        betafactor,
        iotas,
        gammah,
        muh,
        musratio,
        numexclusivepars
    };

    static std::vector<hpartial::params> exclusivepars()
    {
           return { hpartial::betafactor,
               hpartial::iotas,
               hpartial::gammah,
               hpartial::muh,
               hpartial::musratio
           };
    }

    hcohortseir(unsigned numcohorts, const vector<double>& pop) :
        cohortseir<hpartial>(numcohorts,pop)
    {
    }

    virtual std::vector<std::vector<double>> pars2pars(
            unsigned t, const vector<double>& params, const G& g) const
    {
        static std::vector<hpartial::params> computedpars
               = { hpartial::alpha,
                   hpartial::pi,
                   hpartial::eta,
                   hpartial::theta,
                   hpartial::prebeta };

        assert(hpartial::numparams==
                                computedpars.size()
                               +numcommonandcompparams-numcomputationparams
                               +numexclusivepars);


        assert(params.size()==  numcommonandcompparams
                              + numcohorts() * numexclusivepars);

        std::vector<std::vector<double>>
           res(numcohorts());

        vector<double> preparams(hpartial::numparams,0);

        preparams[hpartial::pi] = min(1.0,params[pi] * g.Z(t,PDET));

        double shiftedt = max(t-params[paqshift],0.0);
        unsigned paqtl = static_cast<unsigned>(shiftedt);
        unsigned paqth = static_cast<unsigned>(shiftedt+1);
        double lweight = paqth - shiftedt;

        preparams[hpartial::eta] =
                params[eta0]
//                *  g.Z(t,DAYADJUST)
                ;
        preparams[hpartial::theta] =
                (params[theta0] + params[thetacoef] * g.Z(t,PDET))
//                * g.Z(t,DAYADJUST)
                ;

        double ba =lweight * g.Z(paqtl,BRIGITATTACK) + (1-lweight)* g.Z(paqth,BRIGITATTACK);
        double ba21 = paqtl < 22 ? 0 : (lweight * g.Z(paqtl-14,BRIGITATTACK) + (1-lweight)* g.Z(paqth-14,BRIGITATTACK));

        double c = 1.0 + params[cycleeffect]*cos((g.abstime(shiftedt)-325)/365 * 2*3.141592653);
        double be = 1.0 + params[brigiteeffect] * ba;
        preparams[hpartial::prebeta]
           = c * be * (lweight*g.Z(paqtl,REDUCTIONMEAN) + (1-lweight)*g.Z(paqth,REDUCTIONMEAN))
//              * (lweight*g.Z(paqtl,BFACTOR) + (1-lweight)*g.Z(paqth,BFACTOR))
                * exp(
                    -params[omega] *(lweight*g.Z(paqtl,FEAR) + (1-lweight)*g.Z(paqth,FEAR))
                    -params[omega2] *(lweight*g.Z(paqtl,BETAFACTOR) + (1-lweight)*g.Z(paqth,BETAFACTOR))
                    );


        auto cp = commonpars();
        unsigned srcc=numcomputationparams;
        for(unsigned i=0; i<cp.size(); i++)
            preparams[cp[i]]=params[srcc++];

        auto ep = exclusivepars();
        for(unsigned c=0; c<numcohorts(); c++)
        {
            res[c] = preparams;
            for(unsigned i=0; i<ep.size(); i++)
                res[c][ep[i]] = params[srcc++];
            res[c][hpartial::alpha] =
0.179; // tbd
            res[c][hpartial::iotas] *= 1+ ba21 * params[c<2 ? hybrigitefect : hobrigitefect];
            double pf = 0, ps = 0;
            if(t > 20)
            {
                if(numcohorts()==1)
                {
                    pf = g.Z(t-14,FALL);
                    ps = g.Z(t-7,SALL);
                }
                else
                {
                    pf = g.Z(t-14,F0 + c);
                    ps = g.Z(t-7,S0 + c);
                }
            }
            double vred = ps * 0.1 + (pf-ps) * 0.3 + (1 - pf) * 1;
//clog << t << ": c= " << c << " ps=" << ps << " vred=" << vred << endl;
            res[c][hpartial::prebeta] *= vred;
        }
        assert(srcc==params.size());
        return res;
    }
    virtual double vb(const vector<double>& params ) const
    {
        return params[varbfactor];
    }

    virtual double vp(const vector<double>& params, unsigned i) const
   {
      static int p[hpartial::numstates ] =
    //    E,    Ia,   Ip,   Is,   Iu, R,
       {  ce,   ce,   ce,   ce,   cu, -1,
    //    Edelta, Iadelta, Ipdelta,  Isdelta,  Rdelta,
          cd,     cd,      cd,       cd,       -1,

     //  Hdelta,    Rhdelta,   Ddelta,    Dhdelta,
         ch, -1,        -1,        -1,

     //   Isd, Rd,      Hd,        Rhd,       Dd,     Dhd,
          cd,  -1,      ch,        -1,        -1,     -1 };
       unsigned ci = i % partial().k();

       if(p[ci]==-1)
           return numeric_limits<double>::infinity();
       else
           return params[p[ci]];
    }

};

template <bool single>
class datareader
{
public:
    enum obscolumns {R0,	R20,	R65,	R80,
                     D0,	D20,	D65,	D80,
                     H0,	H20,    H65,    H80,
                     RH0,	RH20,   RH65,   RH80,
                     DH0,	DH20,   DH65,   DH80,
                     RA,	RS,
                      numobscolumns};

    virtual vector<string> obslabels() const = 0;

    virtual unsigned numobs() const = 0;
    virtual vector<double> obs2obs(const vector<double>&) const = 0;

    seirdata read(csv<','>& c, unsigned lag = 0, csv<','>* vc = 0)
    {
        seirdata res;
        bool rfin = false;

        auto lbls = obslabels();
        for(unsigned j=0; j<numobs(); j++)
            res.ylabels.push_back(lbls[j]);

        for(unsigned j=0; j<hcohortseir::numexcolumns; j++)
            res.zlabels.push_back(c(0,numobscolumns+j+1));

//        vector<double> lastobs(numobscolumns,0);

        for(unsigned i=lag; i< c.r()-1; i++)
        {
            res.dates.push_back(c(i+1,0));
            if(c(i+1,1) != "" && !rfin)
            {
                vector<double> origobs;
                for(unsigned j=0; j<numobscolumns; j++)
                    origobs.push_back( c.getdouble(i+1,1+j));
                auto obs = obs2obs(origobs);
                assert(obs.size()==numobs());
                res.y.push_back(dv(obs));
                // lastobs = origobs;
            }
            else
                rfin = true;
            vector<double> exp;
            for(unsigned j=0; j<hcohortseir::numexcolumns; j++)
                exp.push_back(c.getdouble(i+1,1+numobscolumns+j));
            res.z.push_back(dv(exp));            

            vector<double> vac;
            if(single)
                vac.push_back(vc->getdouble(i+1,5));
            else
                for(unsigned j=0; j<4; j++)
                    vac.push_back(vc ? vc->getdouble(i+1,j+1) : 0);
            res.v.push_back(dv(vac));
        }
        res.lag = lag;
        return res;
    }
};


template <typename S, estimationmethod method=emwls>
class hestseir : public S, public estimableseirfilter<method>
{
public:
    hestseir(const uncertain& ax0 = uncertain(dvector(0))) : S(),
        estimableseirfilter<method>(),
        fx0(ax0)
    {}
    virtual uncertain X0(const vector<double>& params) const
    {
        if(fx0.x().size()==0)
            return seirfilter::X0(params);
        else
        {
            assert(fx0.x().size() == this->k());
            return fx0;

/*
            double is= params[initscale];
            auto x= fx0.x();
            auto V = fx0.var();
            auto a = hpartial::numactives;
            for(unsigned i=0; i<a; i++)
                x[i] *= is;
            return { x, block(V.block(0,0,a,a)*is*is,
                              V.block(0,a,a,k()-a) * is,
                              V.block(a,0,k()-a,a) * is,
                              V.block(a,a,k()-a,k()-a)) }; */

        }
    }

private:
    uncertain fx0;
};




#endif // HSEIR_HPP
