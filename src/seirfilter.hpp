#ifndef SEIRFILTER_HPP
#define SEIRFILTER_HPP


#include "orpp.hpp"
#include "orpp/matrices.hpp"
#include "orpp/gnuplot.hpp"
#include "orpp/estimation.hpp"
#include <ostream>

using namespace orpp;
using namespace std;

struct seirdata
{

    vector<string> ylabels;
    vector<string> zlabels;

    vector<string> dates;
    vector<dvector> y;
    vector<dvector> z;
    vector<dvector> v;

    unsigned lag = 0;
    unsigned abstime(unsigned t) const { return t+lag; }

    void output(ostream& os)
    {
        os <<"date";
        for(unsigned i=0; i<ylabels.size(); i++)
            os << "," << ylabels[i];
        for(unsigned i=0; i<zlabels.size(); i++)
            os << "," << zlabels[i];
        os << endl;
        for(unsigned j=0; j<dates.size(); j++)
        {
            os << dates[j];
            for(unsigned i=0; i<ylabels.size(); i++)
            {
                os << "," ;
                if(j<y.size())
                    os << y[j][i];
            }
            for(unsigned i=0; i<zlabels.size(); i++)
            {
                os << ",";
                if(j<z.size())
                    os << z[j][i];
            }
            os << endl;
            // tbd v
        }
    }

    void removebeginning(unsigned toremove)
    {
        dates.erase(dates.begin(), dates.begin() + toremove);
        y.erase(y.begin(), y.begin() + min(toremove,static_cast<unsigned>(y.size())));
        v.erase(v.begin(), v.begin() + min(toremove,static_cast<unsigned>(v.size())));
        z.erase(z.begin(), z.begin() + min(toremove,static_cast<unsigned>(z.size())));
    }
};

class seirfilter
{
    static bool isVar(const dmatrix& W, const string& label, unsigned t,
                      dvector params)
    {
return true;
        for(unsigned j=0; j<W.cols();j++)
        {
            if(W(j,j) < -1e-5)
            {
                cerr << label << "("<<j<<","<<j<<")=" << W(j,j) << endl;
                cerr << label << endl;
                cerr << m2csv(W) << endl;
                cerr << "t=" << t << endl;
                cerr << "pars" << endl;
                cerr << params;
                return false;
            }
            for(unsigned i=j+1; i<W.rows(); i++)
            {
                if(fabs(W(i,j)-W(j,i))>0.003)
                {
                    cerr << label << "("<<j<<","<<i<<")=" << W(i,j)
                         << " != " << W(j,i) << endl;
                    cerr << label << endl;
                    cerr << W << endl;
                    cerr << "t=" << t << endl;
                    cerr << "pars" << endl;
                    cerr << params;
                    return false;
                }
                if(0 && W(i,j) > sqrt(W(i,i)*sqrt(W(j,j))))
                {
                    cerr << label << "("<<i<<","<<j<<")=" << W(i,j)
                         << " > " << sqrt(W(i,i)*W(j,j)) << endl;
                    cerr << label << endl;
                    cerr << W << endl;
                    cerr << "t=" << t << endl;
                    cerr << "pars" << endl;
                    cerr << params;
                    return false;
                }
            }
        }
        return true;
    }

public:

    struct G : private seirdata
    {
        G(const seirdata& d, const seirfilter& am, unsigned evalsize,
            unsigned alongpredlag, bool adiflongpreds, unsigned aestoffset) :
            seirdata(d), longpredlag(alongpredlag), diflongpreds(adiflongpreds),
            estoffset(aestoffset),
            pred(evalsize), predlong(evalsize), est(evalsize),
            hatBs(evalsize),Ps(evalsize),
            contrasts(d.y.size(),0),icontrasts(d.y.size()), contrast(0), sm(am)
           { assert(d.z.size()>=d.y.size()); }

        unsigned longpredlag;
        bool diflongpreds;
        unsigned estoffset;
        vector<uncertain> pred;
        vector<uncertain> predlong;

        vector<uncertain> est;
        vector<dmatrix> hatBs;
        vector<dmatrix> Ps;
        vector<double> contrasts;
        vector<vector<double>> icontrasts;
        double contrast;
        const seirfilter& sm;

        string date(unsigned i) const
        {
            return dates[i];
        }

        dmatrix T(unsigned i) const
        {
            assert(i<Ps.size());
            assert(i<hatBs.size());
            return Ps[i]+hatBs[i];
        }

        unsigned n() const
        {
            assert(y.size());
            return y[0].size();
        }

        unsigned p() const
        {
            assert(z.size());
            return z[0].size();
        }
        virtual dvector Y(unsigned t) const
        {
            assert(t<y.size());
            return y[t];
        }

        double Y(unsigned t, unsigned i) const
        {
            return Y(t)[i];
        }

        unsigned Ysize() const
        {
            return y.size();
        }

        virtual dvector Z(unsigned t) const
        {
            return z[t<z.size() ? t : z.size()-1];
        }

        double Z(unsigned t, unsigned i) const
        {
            return Z(t)[i];
        }

        unsigned Zsize() const
        {
            return z.size();
        }

        double V(unsigned t, unsigned i) const
        {
            return v[t<v.size() ? t : v.size()-1][i];
        }

        struct forecastrecord
        {
            double pred = numeric_limits<double>::quiet_NaN();
            double stderr = numeric_limits<double>::quiet_NaN();
            double act = numeric_limits<double>::quiet_NaN();
            string lbl = "";
        };

        using tforecasts = vector<vector<forecastrecord>>;
        tforecasts forecasts(dmatrix additional = dmatrix())
        {
            vector<vector<forecastrecord>> res;
            unsigned m=additional.rows();
            for(unsigned s=0; s<predlong.size(); s++)
            {
                vector<forecastrecord> r(m);
                if(s>= longpredlag)
                {
                    if(s < y.size())
                    {
                        dvector e = additional * y[s];
                        if(diflongpreds)
                            e -= additional * y[s-longpredlag];

                        for(unsigned i=0; i < additional.rows(); i++)
                            r[i].act = e[i];
                    }

                    uncertain pr = predlong[s];
                    dvector p = additional * pr.x().block(sm.k(),0,sm.n(),1);
                    dmatrix v = additional
                            * pr.var().block(sm.k(),sm.k(),sm.n(),sm.n())
                            * additional.transpose();
                    for(unsigned i=0; i < m; i++)
                        if(!isnan(p[i]))
                        {
                            r[i].pred = p[i];
                            r[i].stderr=sqrt(v(i,i));
                            r[i].lbl = dates[s];
                        }
                }

                res.push_back(r);
            }
            return res;
        }

        struct tforecasteval
        {
            double stdmean;
            double stderr;
            double numoutsideerr;
        };



        static vector<tforecasteval> moments(const vector<vector<forecastrecord>>& r)
        {
            vector<tforecasteval> res;
            for(unsigned j=0; j<r[0].size(); j++)
            {
                double s=0;
                double s2=0;
                unsigned n = 0;
                unsigned noutside =0 ;
                for(unsigned i=0; i<r.size(); i++)
                {
                    if(!(isnan(r[i][j].pred) || isnan(r[i][j].act)))
                    {
                        double dif = (r[i][j].act - r[i][j].pred)/r[i][j].stderr;
                        s += dif;
                        s2 += dif * dif;
                        n++;
                        if(fabs(dif) > 1)
                            noutside++;
                    }
                }
                double mean = s / n;
                res.push_back(
                { mean, sqrt(s2 / n - mean*mean), static_cast<double>(noutside) / n}
                            );
            }
            return res;
        }

        static void fcststognuplot(const vector<vector<forecastrecord>>& r, vector<string> labels)
        {
            for(unsigned i=0; i<r[0].size(); i++)
            {
                string glabel = labels[i] + "_for";
                gnuplot g(glabel);
                for(unsigned j=0; j<r.size();j++)
                    g.datfile() << j << " " << r[j][i].act << " "
                      << r[j][i].pred << " " << r[j][i].stderr << endl;

                g.script()
                << "set style data lines" << endl
                << "set title \"" + labels[i] + " forecasts\"" << endl
                << "plot '" << g.datfn() << "'"
                << " using 2 with lines, '' using 1:3:4 with errorlines" << endl;

                g.process();
            }
        }


        void output(ostream& o,
                    dmatrix additional = dmatrix(),
                    vector<string> lbls = vector<string>())
        {
            assert(ylabels.size()==sm.n());
            o << "date,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_pred,";
            for(unsigned i=0; i<sm.n(); i++)
                o << ylabels[i] << "_pred,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_pred_stdev,";
            for(unsigned i=0; i<sm.n(); i++)
                o << ylabels[i] << "_pred_stdev,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_est,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_est_stdev,";
            for(unsigned i=0; i<sm.n(); i++)
                o << ylabels[i] << ",";
            for(unsigned i=0; i<zlabels.size(); i++)
                o << zlabels[i]<<",";

            for(unsigned i=lbls.size(); i<additional.rows(); i++)
            {
                ostringstream ss;
                ss << "A" << i << endl;
                lbls.push_back(ss.str());
            }
            ostringstream ss;
            if(diflongpreds)
                ss << "(" << longpredlag << ")";

            for(unsigned i=0; i < additional.rows(); i++)
                o << lbls[i] << ss.str() << " act,";
            for(unsigned i=0; i < additional.rows(); i++)
                o << lbls[i] << ss.str()<< " pred,";
            for(unsigned i=0; i < additional.rows(); i++)
                o << lbls[i] << ss.str() << " stdev,";
            o << "contrast,";
            for(unsigned i=0; i<sm.n(); i++)
                o << "C_" + ylabels[i] << ",";
            for(unsigned i=0; i<lag; i++)
                o << endl;
            o << endl;
            for(unsigned s=0; s<pred.size(); s++)
            {
                if(s<dates.size())
                    o << dates[s] << ",";
                else
                    o << ",";
                const uncertain& u =  predlong[s];
                for(unsigned i=0; i<u.dim(); i++)
                    o << u.x()[i] << ",";
                for(unsigned i=0; i<u.dim(); i++)
                    o << sqrt(u.var()(i,i)) << ",";

                if(s < est.size())
                {
                    const uncertain& u = est[s];
                    for(unsigned i=0; i<u.dim(); i++)
                        o << u.x()[i] << ",";
                    for(unsigned i=0; i<u.dim(); i++)
                        o << sqrt(u.var()(i,i)) << ",";
                }
                else
                {
                    for(unsigned i=0; i<est[0].dim(); i++)
                        o << ",,";
                }
                if(s < y.size())
                {
                    for(unsigned i=0; i<y[s].size(); i++)
                        o<< y[s][i] << ",";
                }
                else
                    for(unsigned i=0; i<y[0].size(); i++)
                        o << ",";
                for(unsigned i=0; i<z[0].size(); i++)
                     o << Z(s,i) << ",";
                if(s < y.size() && s>= longpredlag)
                {
                    dvector e = additional * y[s];
                    if(diflongpreds)
                        e -= additional * y[s-longpredlag];

                    for(unsigned i=0; i < additional.rows(); i++)
                        o << e[i] << ",";
                }
                else
                    for(unsigned i=0; i < additional.rows(); i++)
                        o << ",";
                uncertain pr = predlong[s];
                dvector p = additional * pr.x().block(sm.k(),0,sm.n(),1);
                dmatrix v = additional
                        * pr.var().block(sm.k(),sm.k(),sm.n(),sm.n())
                        * additional.transpose();
                for(unsigned i=0; i < additional.rows(); i++)
                    o << p[i] << ",";
                for(unsigned i=0; i < additional.rows(); i++)
                    o << sqrt(v(i,i)) << ",";

                if(s<contrasts.size())
                {
                    o << contrasts[s] << ",";
                    for(unsigned i=0; i<icontrasts[s].size();i++)
                        o << icontrasts[s][i] << ",";
                }
                o << endl;
            }
        }

        vector<double> contrastvars()
        {
            vector<double> s(ylabels.size(),0);
            vector<double> s2(ylabels.size(),0);
            vector<unsigned> n(ylabels.size(),0);

            for(auto r: icontrasts)
            {
                for(unsigned i=0; i<r.size();i++)
                {
                    auto x = r[i];
                    s[i] += x;
                    s2[i] += x * x;
                    n[i]++;
                }
            }
            vector<double> res;
            for(unsigned i=0; i<s.size();i++)
            {
                double m= s[i] / n[i];
                res.push_back( s2[i] / n[i] - m*m);
            }
            return res;
        }

        void report(const string& fn)
        {
            ofstream o(sys::outputfolder()+fn+".tex");
            if(!o)
                throw "cannot open " + sys::outputfolder()+fn;

            o << "\\documentclass[12pt,english,a4paper]{article}" << endl
              << "\\usepackage{graphicx}" << endl
              << "\\begin{document}" << endl << endl;

            o << "\\section*{Residuals}" << endl;



            for(unsigned i=0; i<ylabels.size(); i++)
            {

                o << "\\subsection*{" + ylabels[i] << "}" << endl;

                using namespace std;
                string forlabel = fn+ylabels[i] + "_for";
                string reslabel = fn+ylabels[i] + "_res";
                string stdreslabel = fn+ylabels[i] + "_stdres";
                o << "\\begin{tabular}{ccc}" << endl
//                  <<  "\\includegraphics[width=4cm]{" << forlabel << ".eps}"
//                  <<  "& \\includegraphics[width=4cm]{" << reslabel << ".eps}"
                   <<  "& \\includegraphics[width=6cm]{" << stdreslabel << ".eps}"
                 << "\\end{tabular}" << endl << endl << endl;

                gnuplot g(stdreslabel);

                g.script()
                << "set style data lines" << endl
                << "set zeroaxis" << endl
                << "set title \"" + ylabels[i] + " std resituals\"" << endl
                << "plot '" << g.datfn() << "'"
                << " using 2 with lines" << endl;


                for(unsigned s=0; s<y.size(); s++)
                {
                    auto yy = y[s][i];
                    if(diflongpreds && s > longpredlag)
                        yy -= y[s-longpredlag][i];
                    g.datfile() << s << " "
                      << (yy - predlong[s].x()[i+sm.k()]) / sqrt(predlong[s].var()(i+sm.k(),i+sm.k())) << endl;
                }
                g.process();
            }



            o << "\\end{document}" << endl;
        }



        unsigned abstime(unsigned t) const { return seirdata::abstime(t);}
    };



    virtual unsigned k() const = 0;
    virtual unsigned n() const = 0;
    virtual string statelabel(unsigned i) const
    {
        assert(i < k());
        return std::to_string(i);
    }

    virtual dmatrix P(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dmatrix hatB(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dmatrix F(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dvector I(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dvector gamma(const dvector& x, const dvector& y, unsigned t , const vector<double>& params , const G& g) const
    {
        dvector ret= gamma0(t,params,g) + gamma1(t,params,g) * x + gamma2(t,params,g) * y;
        return ret;
    }
    virtual dvector gamma0(unsigned , const vector<double>&, const G&) const
    {
        dvector ret(n());
        ret.setZero();
        return ret;
    }
    virtual dmatrix gamma1(unsigned , const vector<double>&, const G&) const
    {
        dmatrix ret(n(),k());
        ret.setZero();
        return ret;
    }

    virtual dmatrix gamma2(unsigned , const vector<double>&, const G&) const
    {
        dmatrix ret(n(),k());
        ret.setZero();
        return ret;
    }

    struct evalparams
    {
        unsigned firstcomputedcontrasttime=1;
        unsigned longpredlag = 1;
        unsigned estoffset = 0;
        bool longpredtocontrast = false;
        bool longpredvars = false;
        unsigned ds = 0;
        bool diflongpreds = false;
        double additionalcontrastweight = 100000;
        vector<bool> yfilter = vector<bool>();
        vector<double> additionalwlsweights = vector<double>();
    };


    virtual double contrastaddition(const vector<double>& /* params */, const G& /*r */, evalparams /* ep */) const { return 0; }
    virtual double vb(const vector<double>&/* params */) const { return 1; }
    virtual double vp(const vector<double>&/* params */, unsigned i) const { return numeric_limits<double>::infinity(); }

    virtual uncertain X0(const vector<double>& /*params */) const
    {
        dvector ret(k());
        ret.setZero();
        return uncertain(ret);
    }

    virtual double contrast(unsigned /* t */, const vector<double>& /* params */, G& /*g*/, evalparams /* ep */ ) const
    {
        return 0;
    }


    virtual void adjust(unsigned /*t*/, dvector&  /* N */) const
    {}

    dmatrix T(unsigned t, const vector<double>& pars, const G& g) const
    {
        return P(t,pars,g)+hatB(t,pars,g);
    }

    dmatrix Phi(unsigned t, const vector<double>& pars, unsigned s, const G& g) const
    {
        dmatrix P = this->P(t,pars,g).transpose();
        unsigned k = this->k();
        dmatrix ret(k,k);
        ret.setZero();
        for(unsigned i=0; i<k; i++)
            for(unsigned j=0; j<k; j++)
            {
                ret(i,j) = -P(s,i)*P(s,j);
                if(i==j)
                    ret(i,i)+=P(s,i);
            }
        return ret;
    }

private:
    double overdcoef( double c, double x, double y) const
    {
        return c == numeric_limits<double>::infinity() ? 1
                      : ( c * x + y) / (c+1);
    }
public:
    dmatrix Lambda(unsigned t, const vector<double>& pars, const dvector& x, const dvector& y, const G& g) const
    {
        dmatrix ret( this->k(), this->k());
        ret.setZero();
        for(unsigned i=0; i<this->k(); i++)
            ret += overdcoef(vp(pars,i),x[i],y[i]) * Phi(t,pars,i,g);

        dmatrix eb = hatB(t,pars,g);
        double vf = vb(pars);
        for(unsigned s=0; s<this->k(); s++)
            for(unsigned i=0; i<this->k(); i++)
                ret(i,i) += x[s] * eb(i,s) * vf;

        return ret;
    }

//    double rho(unsigned t, const G& er,unsigned m) const
//    {
//        return radius(er.T(t).block(0,0,m,m));
//    }
    double Rcp(unsigned t, const G& er, unsigned m) const
    {
        const dmatrix& K = er.hatBs[t].block(0,0,m,m);
        for(unsigned i=0; i<K.rows(); i++)
        {
            for(unsigned j=1; j<K.cols(); j++)
                assert(K(i,j)==0);
        }
        dmatrix p(m,1);
        p.setZero();
        p(0,0) = 1;

        // create identity dmatrix (to orpp;
        dmatrix E(m,m);
        E.setZero();
        for(unsigned i=0; i< m; i++)
            E(i,i)=1;


        auto r = K.transpose() * ((E - er.Ps[t].block(0,0,m,m).transpose()).inverse())*p;
        double s = 0;
        for(unsigned i=0;i< m; i++)
            s += r(i,0);
        return s;
    }

    double R(unsigned t, const G& er) const
    {
        unsigned inf = er.Ps.size();
        assert(t<inf);

        const dmatrix& J = er.hatBs[t];
        for(unsigned i=0; i<J.rows(); i++)
        {
            for(unsigned j=1; j<J.cols(); j++)
                assert(J(i,j)==0);
        }
        dmatrix p(k(),1);
        p.setZero();
        p(0,0) = 1;
        double res = 0;
        double thistau = HUGE_VAL;
        for(unsigned tau = t; tau < inf || thistau > 0.00001 ; tau++)
        {
            auto lasttau = min(tau,inf-1);
            dmatrix r = er.hatBs[lasttau].transpose() * p;
            thistau = 0;
            for(unsigned i=0; i<k(); i++)
                thistau += r(i,0);
            res += thistau;
            p = er.Ps[lasttau].transpose() * p;
        }
        return res;
    }

    template<typename TG=G>
    TG eval(const vector<double>& pars,
            const seirdata& data, evalparams ep) const
    {
//        assert(!(ep.diflongpreds && ep.longpredtocontrast));
//        assert(!(!ep.longpredvars && ep.longpredtocontrast));
        assert(data.y.size());
        assert(ep.firstcomputedcontrasttime > 0);

        if(ep.longpredtocontrast)
            assert(ep.longpredlag <= ep.firstcomputedcontrasttime);

        TG g(data,*this, data.y.size() + ep.ds, ep.longpredlag, ep.diflongpreds, ep.estoffset);

        unsigned horizon = g.Ysize()-1;
        assert(ep.estoffset <= horizon);
        unsigned esthorizon = horizon - ep.estoffset;
        unsigned evalhorizon = horizon + ep.ds;
        dvector nanvector(k()+n());
        nanvector.setConstant(numeric_limits<double>::quiet_NaN());

        for(unsigned i=0; i<g.predlong.size(); i++)
            g.predlong[i] = nanvector;
//tbd overhead
        uncertain x0 = X0(pars);

        auto x0mean = x0.x();
        auto x0var = x0.var();
        auto f = F(0,pars,g);
        g.pred[0]=uncertain(stackv(x0mean, f *x0mean),
                   block(x0var, x0var*f.transpose(),
                    f*x0var, f*x0var*f.transpose()));
        g.est[0]=x0;
        g.hatBs[0]=this->hatB(0, pars, g);
        g.Ps[0]=this->P(0, pars, g);
        unsigned i=1;
        for(; i<= evalhorizon; i++)
        {
            bool computingcontrast
                    = (i >= ep.firstcomputedcontrasttime)
                      && ( i<=esthorizon);

            dmatrix T = g.T(i-1);
            dvector Xold;
            dmatrix Wold;
            if(i<=esthorizon+1)
            {
                Xold = g.est[i-1].x();
                Wold = g.est[i-1].var();
                assert(isVar(Wold,"W1",i-1,dv(pars)));
            }
            else
            {
                Xold = g.pred[i-1].x();

                adjust(i,Xold);

                Xold.conservativeResize(k());
                Wold = g.pred[i-1].var().block(0,0,k(),k());
                assert(isVar(Wold,"W2",i-1,dv(pars)));
            }

            dvector Xplus(Xold);
            dvector y(k());
            for(unsigned j=0; j<k(); j++)
            {
                if(Xplus[j] < 0)
                    Xplus[j] = 0;
                y[j] = Xold[j] * Xold[j] + Wold(j,j);
            }

            dmatrix F = this->F(i,pars,g);

            dvector Xnew = T * Xold + I(i-1,pars,g);
            dvector Ynew = F * Xnew;
            dmatrix Wnew = T * Wold * T.transpose() + Lambda(i-1,pars,Xplus,y,g);

            dvector Yplus = F * Xplus;
            dvector gm = gamma(Xplus,square(Xplus),i-1,pars, g);
//clog << i << "gm " << gm[0] << "," << gm[1] << endl;
            dmatrix V =  block( Wnew, Wnew*F.transpose(),
                                F*Wnew, F*Wnew*F.transpose() + diag(gm));

            assert(isVar(V,"V2",i-1,dv(pars)));

            g.pred[i] = uncertain(stackv(Xnew, Ynew),V);
            unsigned ilong = i + ep.longpredlag - 1;
            if(ep.longpredlag && ilong <= evalhorizon && ilong >= ep.firstcomputedcontrasttime)
            {
                G tmpg = g;
                dvector xpred = Xnew;
                dmatrix W;
                dmatrix Tt;

                if(ep.longpredvars)
                {
                    if(ep.diflongpreds)
                    {
                        W=dmatrix(k(),k());
                        W.setZero();
                        Tt=dmatrix::Identity(k(),k());
                    }
                    else
                        W = Wnew;
                }
                dvector xpredplus = Xplus;

                for(unsigned j=1;j<ep.longpredlag;j++ )
                {
                    tmpg.est[i+j-1] = { xpred, W };

                    auto thisT = this->T(i+j-1,pars,tmpg);
                    xpred = thisT*xpred + I(i+j-1,pars,tmpg);

                    if(ep.longpredvars)
                    {
                        xpredplus = xpred;
                        dvector y(k());

                        for(unsigned j=0; j<k(); j++)
                        {
                            if(xpredplus[j] < 0)
                                xpredplus[j] = 0;
                            y[j] = xpred[j]*xpred[j] + W(j,j);
                        }

                        W = T * W * T.transpose() + Lambda(i+j-1,pars,xpredplus,y,tmpg);
                        if(ep.diflongpreds)
                            Tt *= thisT;
                    }
                }

                auto Flong = this->F(ilong,pars,tmpg);

                dvector lp;
                if(ep.diflongpreds)
                   lp= stackv(xpred-Xold, Flong*xpred-F*Xold);
                else
                   lp = stackv(xpred, F*xpred);
                if(ep.longpredvars)
                {
                    dvector gmnew = gamma(xpredplus,square(xpredplus), ilong-1,pars, tmpg);

                    if(ep.diflongpreds)
                    {
                        auto trsf = Tt-dmatrix::Identity(k(),k());
                        if(i > esthorizon)
                        {
                            W += trsf * Wnew * trsf.transpose();
                            g.predlong[ilong] = {
                                lp,
                                block( W, W*Flong.transpose(),
                                       Flong*W, Flong*W*Flong.transpose()
                                          + diag(gm+gmnew))
                                };
                        }
                        else
                        {
                            auto FTt = F * Tt;
                            g.predlong[ilong] = {
                                lp,
                                block( W+ trsf * Wnew * trsf.transpose(),
                                       W*Flong.transpose() + trsf * Wnew * FTt.transpose(),
                                       Flong*W + FTt * Wnew * trsf.transpose(),
                                       Flong*W*Flong.transpose() + FTt * Wnew * FTt.transpose()
                                          + diag(gmnew))
                                };

                        }
                    }
                    else
                    {
                        g.predlong[ilong] = {
                            lp,
                            block( W, W*Flong.transpose(),
                                   Flong*W, Flong*W*Flong.transpose()
                                      + diag(gmnew))
                            };
                    }
                }
                else
                    g.predlong[ilong] = lp;
            }

            if(i<=esthorizon)
            {
                dmatrix Vxx = V.block(0,0,k(),k());
                dmatrix Vxy = V.block(0,k(),k(),n());
                dmatrix Vyx = V.block(k(),0,n(),k());
                dmatrix Vyy = V.block(k(),k(),n(),n());

                if(Vyy.isZero())
                {
                    g.est[i] = uncertain(Xnew, Vxx); // ?? divný
//                    assert(!computingcontrast);
                }
                else
                {
                    dmatrix Vinv = pseudoinverse(Vyy); // Vyy.inverse();

                    dvector Xest = Xnew + Vxy * Vinv * (g.Y(i)-Ynew);
                    dmatrix West =Vxx + ((-1) * Vxy) * Vinv * Vyx;
                    assert(isVar(West,"W",i-1,dv(pars)));
                    g.est[i] = uncertain(Xest,West);

                    if(computingcontrast)
                    {
                        double c = contrast(i,pars,g,ep);
                        g.contrasts[i] = c;
                        g.contrast += c;
                    }
                }
            }
            else // i<=t
            {
                auto& p = g.pred[i];
                g.est[i] = uncertain(p.x().block(0,0,k(),1),
                                     p.var().block(0,0,k(),k()) );
            }


            g.hatBs[i] = this->hatB(i,pars, g);
            g.Ps[i] = this->P(i,pars, g);
        }

        g.contrast += ep.additionalcontrastweight * contrastaddition(pars,g,ep);
//throw;
        return g;
    }

    virtual dmatrix activesubmatrix( const dmatrix& T) const
    {
        return T;
    }

    double radius(unsigned t, const vector<double>& pars, const G& g) const
    {
        return ::radius(activesubmatrix(T(t,pars,g)));
    }

};


#endif // SEIRFILTER_HPP
