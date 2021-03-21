#ifndef COMPLEXEST_HPP
#define COMPLEXEST_HPP

#include "hfourseir.hpp"
#include "estimableseirfilter.hpp"
/*
template <typename S>
inline vector<seirparaminit> doest(S& seir,
                             const vector<seirparaminit>& params,
                             const seirdata& ad,
                             seirfilter::evalparams eep,
                             double timest)
{

    uncertain u;
    seir.estimate(params, ad, eep, u, timest);
    vector<seirparaminit> res = params;
    for(unsigned i=0; i<res.size(); i++)
    {
        res[i].initial = u.x()[i];
        clog << "{\"" << params[i].name << "\","
            << res[i].initial << "," << params[i].lower << "," << params[i].upper << ","
            << (params[i].omit ? "true" : "false") << "},";
        clog << "// (" << sqrt(u.var()(i,i)) << ") "
                << endl;
    }
    return res;
}

template <typename S>
virtual double complexest(uncertain& x0,
                const vector<seirparaminit>& params,
                const seirdata& ad,
                unsigned estoffset,
                double acweight,
                uncertain& res,
                double timest = numeric_limits<double>::infinity()
                )
{
//        for(unsigned i=0; i<=p.size(); i++)
//            p.omit = false;
//    seirfilter::evalparams eep = {40,7,estoffset,true, false};

//    auto p1 = doest(p,ad,eep,timest / 3.0);



    hestseir<S,emwlsstd> esv(x0);

    vector<bool> yf(esv.n(),false);
    yf[hfourseir::R0] = yf[hfourseir::R20]
            = yf[hfourseir::R65] = yf[hfourseir::R80] = true;

    seirfilter::evalparams vep = {40,7,estoffset,true, true};
    vep.additionalcontrastweight = acweight;
    vep.yfilter = yf;

    auto p = params;
    for(unsigned i=hcohortseir::firstdisp; i<=hcohortseir::lastvar; i++)
          p[i].omit = true;

    vector<unsigned> firstcommonpars = {hcohortseir::ce, hcohortseir::cu,
        hcohortseir::rcoef, hcohortseir::omega,hcohortseir::pi,
           hcohortseir::thetacoef, hcohortseir::eta0,
           hcohortseir::theta0,
           hcohortseir::sigma, hcohortseir::varsigma,
           hcohortseir::ufactor, hcohortseir::vfactor
    };

    for(unsigned i=0; i<firstcommonpars.size(); i++)
        p[i].omit = params[i].omit;

    auto p1 = doest<es>(esv,p,ad,eep,timest / 5.0);

    return 0;
}
*/

#endif // COMPLEXEST_HPP
