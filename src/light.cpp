#include "lseir.hpp"

vector<seirparaminit> cz=
{
//    v, ce, cs, ccoef, dcoef,
//                      alpha, sigma, varsigma, gammas, gammaa, mu, beta, tshift, ralpha, numparams};
    {"v",26,0,2000,true},// (0.431221)
    {"ce",1e+07,0,1e+09,true},// (1.66768e+06)
    {"cs",1e+07,0,1e+09,true},// (2.07871e+08)
    {"ccpef",0,0,100,true},// (1.1647e-06)
    {"dcoef",0,0,100,true},// (9.68045e-05)
    {"alpha",0.179,0,1,true},// (0.0137321)
    {"sigma",0.1787,0,1,true},// (0.00039807)
    {"varsigma",0.2212,0,1,true},// (0.000543095)
    {"gamma",0.296364,0,1,false},// (0.000563102)
    {"mu",0.00572859,1e-10,1,false},// (1.09572e-05)
    {"beta",0.0821299,0,6,false},// (0.00118148)
    {"tshift",0,0,14,true},// (0.000320735)
    {"talpha",0.01,0,0.1,true},// (5.52741e-05)
};

constexpr unsigned numobsforest= 300;
constexpr unsigned fineobsforest= 40;
constexpr unsigned obsomit = 50;


template <estimationmethod M>
double finddisp(estlseir<M>& es, const seirdata& d,
                        const vector<double>& ap,
                        bool s)
{
    dmatrix b(1,2);
    b.setZero();
    b(0,s) = 1;

    unsigned disp = s ? lseir::cs : lseir::ce;

    seirfilter::evalparams oep = {obsomit,7,0,true, true, 0, true};

    double hi= 1e5;
    double lo= 1e-5;
    auto p = ap;
    for(;;)
    {
       p[disp] = hi;
       seirfilter::G g = es.eval(p, d, oep);
       auto ms = seirfilter::G::moments(g.forecasts(b));
       if(ms[0].stderr < 1)
           hi *= 2;
       else
           break;
    }

    for(;;)
    {
       p[disp] = lo;
       seirfilter::G g = es.eval(p, d, oep);
       auto ms = seirfilter::G::moments(g.forecasts(b));
       if(ms[0].stderr > 1)
           lo /= 2;
       else
           break;
    }

    for(;;)
    {
        double mid = (hi + lo) / 2;
        p[disp] = mid;
        seirfilter::G g = es.eval(p, d, oep);
        auto ms = seirfilter::G::moments(g.forecasts(b));
        if(ms[0].stderr > 1)
            hi = mid;
        else
            lo = mid;
        if(hi-lo < 0.01)
            return mid;
    }
};

enum estmode { noest, lastadjust, disperse, initial, full };

template <estimationmethod M>
vector<double> estimate(estlseir<M>& es,  const seirdata& d,
                        const vector<seirparaminit>& pinit,
                        double atime, unsigned numomit)
{

    seirfilter::evalparams eep = {numomit,7,0,
                                  true, M==emwls ? false : true,0, false, 50000};


    uncertain res;

    clog << "Estimating by wls." << endl; // tbd not always

    es.estimate(pinit,d, eep , res, atime);

    auto rp = vd(res.x());

    clog << "{"<< endl;
    for(unsigned i=0; i<pinit.size(); i++)
    {
        clog << "{\"" << pinit[i].name << "\","
            << rp[i] << "," << pinit[i].lower << "," << pinit[i].upper << ","
            << (pinit[i].omit ? "true" : "false") << "},";
        clog << "// (" << sqrt(res.var()(i,i)) << ") "
                << endl;
    }
    clog << "}," << endl;
    return rp;
}


void runlite(const string& country, const string& cabbr,
         const vector<seirparaminit>& apinit, double atime = 3*60, estmode mode = full )
{ // ucitele: ucinovaclock.csv ucihhlllock.csv

    assert(pinit.size() == lseir::numparams);

    ldatareader r;

    seirdata fulld = r.read(
                "/home/martin/tmp/integrace/" + cabbr + ".csv",
                country,
                "/home/martin/Documents/s/covid/data/COVID-19/"
                );
//    fulld.output(clog);

    dvector x0(lseir::numstates);
    x0.setZero();
    x0[lseir::E] = 3;

    estlseir<emwls> es(x0);

    vector<double> rp = spi2v(apinit);
    if(mode == full)
    {
        clog << "Estimating full period" << endl;
        rp = estimate<emwls>(es,fulld,apinit,atime / 3,obsomit);
    }


    clog << "Setting initial values." << endl;


    seirfilter::evalparams iep = {1000,0,0,true, false};

    seirfilter::G rinit = es.eval(rp, fulld, iep);

    unsigned sechor = fulld.y.size()-numobsforest;
    auto d=fulld;
    d.removebeginning(sechor);

    es.fx0 = rinit.est[sechor];

    auto pinit = apinit;
    v2spi(rp,pinit);

    if(mode >= initial)
        rp = estimate<emwls>(es,d,pinit,atime / 3,obsomit);

    dmatrix b(2,2);
    b.setIdentity();



    if(0 && mode >= disperse)
    {
        clog << "setting dispersions" << endl;
        double ce =  finddisp(es,d,rp,false);
        clog << "ce set to " << ce<< endl;
        rp[lseir::ce]=ce;

//        rp = estimate<emwls>(es,d,pinit,atime / 3);

        double cs =  finddisp(es,d,rp,true);
        clog << "cs set to " << cs << endl;

        rp[lseir::cs] = cs;

        v2spi(rp,pinit);
        for(unsigned i=0; i<pinit.size(); i++)
            pinit[i].omit = true;
pinit[lseir::dcoef].omit = false;
pinit[lseir::ccoef].omit = false;
pinit[lseir::ce].omit = false;
pinit[lseir::cs].omit = false;

estlseir<emwlsstd> ses(x0);

        rp = estimate<emwlsstd>(ses,d,pinit, atime / 3, obsomit);
    }

    if(mode >= lastadjust)
    {
        clog << "last adjustment" << endl;
        v2spi(rp,pinit);
        for(unsigned i=0; i<pinit.size(); i++)
            pinit[i].omit = true;
        pinit[lseir::beta].omit=pinit[lseir::gammas].omit=pinit[lseir::mu].omit = false;
pinit[lseir::tshift].omit=false;
        rp = estimate<emwls>(es,d,pinit, atime / 3, d.y.size()
-40
                             );

    }

    seirfilter::evalparams evalep =
      {
d.y.size() - 40
        ,7,0,true, true,ldatareader::forecastlen, true};

    seirfilter::G g = es.eval(rp, d, evalep);
    auto f = g.forecasts(b);
    auto ms = seirfilter::G::moments(f);
    clog << "C: " << ms[0].stderr  << "(" << ms[0].stdmean << ","
                          << ms[0].numoutsideerr << ")" << endl;
    clog << "D: " << ms[1].stderr  << "(" << ms[1].stdmean << ","
                          << ms[1].numoutsideerr << ")" << endl;

    ofstream o(sys::outputfolder()+country + "_outputest.csv");
    vector<string> lbls = {"C","D"};
    g.output(o,b,lbls);

    for(unsigned i=0; i<f.size(); i++ )
        for(unsigned j=0; j<2; j++)
            f[i][j].stderr *= ms[j].stderr;
    seirfilter::G::fcststognuplot(f, {cabbr + "_C", cabbr+"_D"});
};


int main()
{
    try {
        sys::setoutputfolder("../output/");
        runlite("Czechia","CZ",cz);
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


//    clog << "time: " << time(0)-startt << " seconds." << endl;



    return 0;

}
