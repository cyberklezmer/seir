#define GOOGLE

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
#ifdef GOOGLE
    {"retail",0,-1,1,false},
    {"grocery",0,-1,1,true},
    {"parks",0,-1,1,true},
    {"ransit",0,-1,1,false},
    {"workplaces",0,-1,1,false},
    {"residential",0,-1,1,true},
#endif
    {"betafactor",1,0.5,5,true},
    {"mufactor",1,0.5,5,true}
};

constexpr unsigned numobsforest= 300;
constexpr unsigned obsomit = 50;


template <estimationmethod M>
double finddisp(estlseir<M>& es, const seirdata& d,
                        const vector<double>& ap,
                        bool s, unsigned estoffset)
{
    dmatrix b(1,2);
    b.setZero();
    b(0,s) = 1;

    unsigned disp = s ? lseir::cs : lseir::ce;

    seirfilter::evalparams oep = {obsomit,7,estoffset,true, true, 0, true};

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

enum estmode { noest, checkdata, lastadjust, disperse, initial, full };

template <estimationmethod M>
vector<double> estimate(estlseir<M>& es,  const seirdata& d,
                        const vector<seirparaminit>& pinit,
                        double atime, unsigned numomit, unsigned estoffset)
{

    seirfilter::evalparams eep = {numomit,7,estoffset,
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


seirfilter::G::tforecasts runlite(const string& country, const string& cabbr,
         const vector<seirparaminit>& apinit, double atime, estmode mode, unsigned estoffset )
{ // ucitele: ucinovaclock.csv ucihhlllock.csv

    assert(apinit.size() == lseir::numparams);

    ldatareader r;

    seirdata fulld = r.read(
                "../input/country_data/"
                 + cabbr + ".csv",
                country,
                "../input/jhu/"
                );
    for(unsigned i=fulld.y.size()-estoffset; i<fulld.z.size(); i++)
        fulld.z[i] = fulld.z[fulld.y.size()-estoffset-1];
//    fulld.output(clog);
    if(mode == checkdata)
        return seirfilter::G::tforecasts(0);

    dvector x0(lseir::numstates);
    x0.setZero();
    x0[lseir::E] = 3;

    estlseir<emwls> es(x0);

    vector<double> rp = spi2v(apinit);
    if(mode == full)
    {
        clog << "Estimating full period" << endl;
        rp = estimate<emwls>(es,fulld,apinit,atime / 3,obsomit, estoffset);
    }


    clog << "Setting initial values." << endl;


    seirfilter::evalparams iep = {1000,0,estoffset,true, false};

    seirfilter::G rinit = es.eval(rp, fulld, iep);

    unsigned sechor = fulld.y.size()-numobsforest;
    auto d=fulld;
    d.removebeginning(sechor);

    es.fx0 = rinit.est[sechor];

    auto pinit = apinit;
    v2spi(rp,pinit);

    if(mode >= initial)
        rp = estimate<emwls>(es,d,pinit,atime / 3,obsomit,estoffset);

    dmatrix b(2,2);
    b.setIdentity();

    if(0 && mode >= disperse)
    {
        clog << "setting dispersions" << endl;
        double ce =  finddisp(es,d,rp,false,estoffset);
        clog << "ce set to " << ce<< endl;
        rp[lseir::ce]=ce;

//        rp = estimate<emwls>(es,d,pinit,atime / 3, estoffset);

        double cs =  finddisp(es,d,rp,true, estoffset);
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

        rp = estimate<emwlsstd>(ses,d,pinit, atime / 3, obsomit, estoffset);
    }

    if(mode >= lastadjust)
    {
        clog << "last adjustment" << endl;
        v2spi(rp,pinit);
        for(unsigned i=0; i<pinit.size(); i++)
            pinit[i].omit = true;
        pinit[lseir::betafactor].omit=false;
        pinit[lseir::mufactor].omit=false;
        rp = estimate<emwls>(es,d,pinit, atime / 6, d.y.size()-lseir::fineobs-estoffset, estoffset);
    }

    seirfilter::evalparams evalep =
      { static_cast<unsigned>(d.y.size())-lseir::fineobs-estoffset
        ,7,estoffset,true, true,ldatareader::forecastlen, true};


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

    f.resize(f.size()-estoffset);
    f.erase(f.begin(),f.begin()+d.y.size()-estoffset-lseir::fineobs);
    for(unsigned i=0; i<f.size(); i++ )
        for(unsigned j=0; j<2; j++)
            f[i][j].stderr *= ms[j].stderr;
    ostringstream eo;
    eo << estoffset;
    for(unsigned i=0; i<f.size()-ldatareader::forecastlen; i++)
        for(unsigned j=0; j<2; j++)
            f[i][j].pred =f[i][j].stderr = numeric_limits<double>::quiet_NaN();
    seirfilter::G::fcststognuplot(f, {cabbr + "_C"+eo.str(), cabbr+"_D"+eo.str()});
    return f;
};

struct countryinfo
{
    string abbr;
    string name;
};

vector<countryinfo> countries =
{
    {"AT", "Austria"},
    {"BE", "Belgium"},
    {"BG", "Bulgaria"},
    {"CH", "Switzerland"},
    {"CZ", "Czechia"},
    {"DE", "Germany"},
    {"DK", "Denmark"},
    {"EE", "Estonia"},
    {"ES", "Spain"},
    {"FI", "Finland"},
    {"FR", "France"},
    {"GB", "United Kingdom"},
    {"GR", "Greece"},
    {"HR", "Croatia"},
    {"HU", "Hungary"},
    {"IE", "Ireland"},
    {"IT", "Italy"},
    {"LT", "Lithuania"},
    {"LV", "Latvia"},
    {"MT", "Malta"},
    {"NL", "Netherlands"},
    {"NO", "Norway"},
    {"PL", "Poland"},
    {"PT", "Portugal"},
    {"RO", "Romania"},
    {"SE", "Sweden"},
    {"SI", "Slovenia"},
    {"SK", "Slovakia"}
};

void light(bool expforecast)
{
    try {
        sys::setoutputfolder("../output/");

        if(expforecast) // evaluate
        {
            string submdate = "2021-04-05";
            ofstream ecvs(sys::outputfolder()
                          +submdate+"-bisop-seirfilterlite.csv");
            if(!ecvs)
                throw "cannot open forecast file";
            ecvs << "forecast_date,scenario_id,target,target_end_date,location,type,quantile,value" << endl;


//            ofstream o(sys::outputfolder()+"inc.tex");
//            if(!o)
//                throw "cannot open inctex";
            for(auto r : countries)
            {
                try
                {
                    static vector<vector<double>> cfscale =
                    {
                        {0.85, 1.15, 1.8, 2.5},
                        {1.5, 1.7, 1.9, 2.20 }
                    };
                    auto f=runlite(r.name, r.abbr,cz,5*60,full, 0);
                    for(unsigned i=0; i<2; i++)
                    {
                        unsigned t=0;
                        for(;t<f.size();t++)
                        {
                            if(f[t][i].lbl == submdate)
                                break;
                        }
                        if(t==f.size())
                            throw submdate + " not found.";
                        seirfilter::G::tforecasts sf;
                        unsigned k=0;
                        for(int s = t-2 - 4*7; s < t + 4*7; s+=7)
                        {
                            seirfilter::G::forecastrecord fr;
                            fr.act = f[s][i].act;
                            fr.lbl = f[s][i].lbl;
                            if(s>t)
                            {
                                fr.pred = f[s][i].pred;
                                fr.stderr = f[s][i].stderr * 1.96 * cfscale[i][k];

                                static vector<double> qs = {0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.975, 0.99};
                                static vector<double> qsv = {-2.32634787404084,	-1.95996398454005,	-1.64485362695147,	-1.2815515655446,	-1.03643338949379,	-0.841621233572914,	-0.674489750196082,	-0.524400512708041,	-0.385320466407567,	-0.2533471031358,	-0.125661346855074,	0,	0.125661346855074,	0.2533471031358,	0.385320466407568,	0.524400512708041,	0.674489750196082,	0.841621233572914,	1.03643338949379,	1.2815515655446,	1.64485362695147,	1.95996398454005,	2.32634787404084};

                                for(unsigned j=0; j<qs.size()+1; j++)
                                {
                                    ecvs << submdate << ",forecast,"
                                         << k+1 << " wk ahead inc " << (i==0 ? "case," : "death," )
                                         <<  fr.lbl << "," << r.abbr << ",";
                                    if(j<qs.size())
                                        ecvs << "quantile," << qs[j] << "," << round(fr.pred + f[s][i].stderr * cfscale[i][k] * qsv[j]);
                                    else
                                        ecvs << "point,NA," << fr.pred;
                                    ecvs << endl;
                                }
                                k++;
                            }
                            sf.push_back({fr});
                        }

                        static vector<string> lbls={"C", "D"};
                        seirfilter::G::fcststognuplot(sf, {r.abbr + "_" + lbls[i]});
                    }
//                    o << "\\begin{tabular}{ccc}" << endl
//                       <<  r.abbr << "& \\includegraphics[width=4cm]{" << r.abbr << "_C_for.eps}"
//                       << "& \\includegraphics[width=4cm]{" << r.abbr << "_D_for.eps}"
//                       << "\\end{tabular}" << endl;


                }
                catch (std::exception& e)
                    {
//                        o << r.abbr << " " << e.what() << endl;
                    }
                    catch (const char* m)
                    {
//                        o << r.abbr << " " << m << " " << endl;
                    }
                    catch (const string& m)
                    {
//                        o << r.abbr << " " << m << " " << endl;
                    }
                    catch(...)
                    {
//                        o << r.abbr << " unknown exception" << endl;
                    }
//               o << endl;
            }
        }
        else
        {
            ofstream otex(sys::outputfolder()+"eval.tex");
            if(!otex)
                throw "cannot open otex";


            ofstream o(sys::outputfolder()+"eval.csv");
            if(!o)
                throw "cannot open evalcsv";
            for(auto r : countries)
            {
                try
                {
                    otex << "\\begin{tabular}{ccc}" << endl;

                    o << r.abbr;
                    vector<unsigned> cfails(5,0);
                    vector<unsigned> dfails(5,0);
                    unsigned trials = 0;
                    for(unsigned t=35;
trials < 5 &&
                        t < numobsforest / 2; t += 14, trials++)
                    {
                        auto f = runlite(r.name, r.abbr,cz,
20,
 //                                       5*60,
                                         full, t);
                        for(unsigned i=0;i<=5;i++)
                        {
                            auto rc=f[f.size()-1-i*7][0];
                            auto rd=f[f.size()-1-i*7][1];
                            if(!(isnan(rc.pred) || isnan(rc.act) || isnan(rc.stderr)))
                            {
                                if(fabs(rc.act-rc.pred) > rc.stderr)
                                    cfails[i]++;
                            }
                            if(!(isnan(rd.pred) || isnan(rd.act) || isnan(rd.stderr)))
                            {
                                if(fabs(rd.act-rd.pred) > rd.stderr)
                                    dfails[i]++;
                            }
                        }

                        otex <<  r.abbr <<  " " << t
                              << "& \\includegraphics[width=3cm]{" << r.abbr << "_C"
                              << t << "_for.eps}"
                              << "& \\includegraphics[width=3cm]{" << r.abbr << "_D"
                              << t << "_for.eps} \\\\";
                    }
                    for(unsigned i=0; i<5; i++)
                        o << "," << cfails[i];
                    for(unsigned i=0; i<5; i++)
                        o << "," << dfails[i];
                    o <<  "," << trials << endl;

                    otex << "\\end{tabular}" << endl;

                }
                catch (std::exception& e)
                {
                    o << r.abbr << " " << e.what() << endl;
                }
                catch (const char* m)
                {
                    o << r.abbr << " " << m << " " << endl;
                }
                catch (const string& m)
                {
                    o << r.abbr << " " << m << " " << endl;
                }
                catch(...)
                {
                    o << r.abbr << " unknown exception" << endl;
                }
            }

        }
    }
    catch (std::exception& e) {
        std::cerr << e.what() << endl;
    }
    catch (const char* m)
    {
           cerr << m << endl;
    }
    catch (const string& m)
    {
           cerr << m << endl;
    }
    catch(...)
    {
        cerr << "unknown exception" << endl;
    }


//    clog << "time: " << time(0)-startt << " seconds." << endl;
}
