#ifndef SEIR_HPP
#define SEIR_HPP

#include "orpp.hpp"
#include "orpp/matrices.hpp"
#include "orpp/estimation.hpp"

using namespace orpp;
using namespace std;

inline bool isVar(const dmatrix& W, const string& label, unsigned t,
                  dvector params)
{
    for(unsigned j=0; j<W.cols();j++)
    {
        if(W(j,j) < -1e-5)
        {
            cerr << label << "("<<j<<","<<j<<")=" << W(j,j) << endl;
            cerr << label << endl;
            cerr << W << endl;
            cerr << "t=" << t << endl;
            cerr << "pars" << endl;
            cerr << params;
            return false;
        }
        for(unsigned i=j+1; i<W.rows(); i++)
        {
            if(fabs(W(i,j)-W(j,i))>0.001)
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

class seirmodel
{
public:
    struct history { vector<dvector> y; vector<dvector> z; };
protected:
    struct G: public history
    {
        vector<dvector> est;
    };

public:

    enum contrasttype {ectll, ectwls};
    struct data: public history
    {
        vector<string> dates;
        vector<string> ylabels;
        vector<string> zlabels;
        vector<vector<double>> w;
    };
protected:
    dmatrix T(unsigned t, const vector<double>& pars, const G& g) const
    {
        return P(t,pars,g).transpose()+J(t,pars,g).transpose();
    }
    dmatrix Phi(unsigned t, const vector<double>& pars, unsigned s, const G& g) const
    {
        dmatrix P = this->P(t,pars,g);
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

    dmatrix bigPhi(unsigned t, const vector<double>& pars, unsigned s,const G& g) const
    {
        dmatrix phi = Phi(t,pars,s,g);
        return block(phi,phi*H.transpose(),H*phi,H*phi*H.transpose());
    }


    dmatrix bigIwGamma(unsigned t, const vector<double>& pars,
                const dvector& N, const G& g ) const
    {
        dvector n(N.size());
        for(unsigned i=0; i<N.size(); i++)
            n[i] = max(N[i],0.0);

        dmatrix ret(k()+m(),k()+m());
        ret.setZero();
        dvector Jx = J(t,pars,g).transpose() * n;
        for(unsigned i=0; i<Jx.size(); i++)
            if(Jx[i]<0.0)
                throw "negative J";
        dmatrix D= varfactor(pars) * diag(Jx);
dvector nn(n);
for(unsigned i=0; i<n.size(); i++)
    nn[i] *= nn[i];

        return block( D, D*H.transpose(), H*D, H*D*H.transpose()
                      + diag(Gamma(t,pars, g)*
nn // n
                             ));
    }

    virtual double contrast(contrasttype act, const dvector& o,
                            const dvector& My,
                            const dmatrix& Vyy, const dmatrix& Vinv, vector<double> w) const
    {
        if(act == ectll)
            return -(o.size()/2.0) * ::log(2*3.1415926)
            -0.5 * ::log(Vyy.determinant())
            -0.5 * ((o-My).transpose()*Vinv*(o-My))(0,0);
        else if(act == ectwls)
        {
            double res = 0;
            for(unsigned i=0; i<My.size(); i++)
            {
                double x = o[i]-My[i];
                res += w[i] * x * x;
            }
            return -res;
        }
        throw "unknown contrast method";
    }


public:

    struct evalresult
    {
        evalresult(unsigned t, unsigned s, const seirmodel& am) :
            pred(s+1), est(t+1),Ts(s+1),Js(s+1),Ps(s+1),
            contrasts(t,0), contrast(0), sm(am)
           { assert(s>=t); }
        vector<uncertain> pred;
        vector<uncertain> est;
        vector<dmatrix> Ts;
        vector<dmatrix> Js;
        vector<dmatrix> Ps;
        vector<double> contrasts;
        double contrast;
        const seirmodel& sm;
        void output(ostream& o, const data& d)
        {
            assert(d.ylabels.size()==sm.m());
            o << "date,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_pred,";
            for(unsigned i=0; i<sm.m(); i++)
                o << d.ylabels[i] << "_pred,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_pred_var,";
            for(unsigned i=0; i<sm.m(); i++)
                o << d.ylabels[i] << "_pred_var,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_est,";
            for(unsigned i=0; i<sm.k(); i++)
                o << sm.statelabel(i)<<"_est_var,";
            for(unsigned i=0; i<sm.m(); i++)
                o << d.ylabels[i] << ",";
            for(unsigned i=0; i<d.zlabels.size(); i++)
                o << d.zlabels[i]<<",";
            o << endl;
            for(unsigned s=0; s<pred.size(); s++)
            {
                o << d.dates[s] << ",";
                const uncertain& u = pred[s];
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
                if(s < d.y.size() && s < d.z.size())
                {
                    for(unsigned i=0; i<d.y[s].size(); i++)
                        o<< d.y[s][i] << ",";
                    for(unsigned i=0; i<d.z[s].size(); i++)
                        o<< d.z[s][i] << ",";
                }
                o << endl;
            }
        }
    };

    double rho(unsigned t, const evalresult& er,unsigned i) const
    {
        return radius(er.Ts[t].block(0,0,i,i));
    }
    double Rcp(unsigned t, const evalresult& er, unsigned m) const
    {
        const dmatrix& K = er.Js[t].block(0,0,m,m);
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

    double R(unsigned t, const evalresult& er) const
    {
//clog << "Computing R from at " << t << endl;
        unsigned inf = er.Ps.size();
        assert(t<inf);

        const dmatrix& J = er.Js[t];
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
            dmatrix r = er.Js[lasttau].transpose() * p;
            thistau = 0;
            for(unsigned i=0; i<k(); i++)
                thistau += r(i,0);
            res += thistau;
//clog << "lv:" << lastval << endl;
//clog <<er.Js[tau] << endl;
//clog << p << endl << endl;
//clog << lastval << endl;
            p = er.Ps[lasttau].transpose() * p;
//clog << "t=" << t << endl;
//clog << r << endl;
        }
        return res;
    }

    virtual void adjust(unsigned /*t*/, dvector&  /* N */) const
    {}

    template <bool contrasts=true>
    evalresult eval(const vector<double>& pars,
                    const data& d, unsigned ds=0, contrasttype act = ectll) const
    {
        assert(d.y.size());
        unsigned t = d.y.size()-1;
        assert(ds >= 0);
        assert(d.z.size()>=t+ds);
        evalresult res(t,t+ds,*this);
        G g;
        g.y = {d.y[0]};
        g.z = {d.z[0]};
        dvector x0 = dv(X0(pars,g));
        res.pred[0]=stackv(x0, H*x0);
        g.est = {x0};
        res.est[0]=x0;
        res.Ts[0]=this->T(0, pars, g);
        res.Js[0]=this->J(0, pars, g);
        res.Ps[0]=this->P(0, pars, g);
        unsigned i=1;
        for(; i<= t+ds; i++)
        {
            dmatrix& T = res.Ts[i-1];

            dmatrix THT = stackm(T,H*T);

            dvector N;
            dmatrix W;
            if(i<=t+1)
            {
                N = res.est[i-1].x();
                W = res.est[i-1].var();
                assert(isVar(W,"W1",i-1,dv(pars)));
            }
            else
            {
                N = res.pred[i-1].x();

                adjust(i,N);

                N.conservativeResize(k());
                W = res.pred[i-1].var().block(0,0,k(),k());
                assert(isVar(W,"W2",i-1,dv(pars)));
            }
            dvector Mx = T * N + dv(I(i-1,pars,g));
            dvector My = H * Mx;

            dmatrix V = THT * W * THT.transpose();

//cerr << "W" << endl << W << endl;
//cerr << "THT" << endl << THT << endl;
//cerr << "THTT" << endl << THT.transpose() << endl;

            assert(isVar(V,"V2",i-1,dv(pars)));
            for(unsigned j=0; j<k(); j++)
                V = V + max(N[j],0.0) * bigPhi(i-1,pars,j,g);
            assert(isVar(V,"V3",i-1,dv(pars)));

            V = V + bigIwGamma(i-1, pars, N,g);
            assert(isVar(V,"V4",i-1,dv(pars)));

            res.pred[i] = uncertain(stackv(Mx, My),V);
            if(i<=t)
            {
                dmatrix Vxx = V.block(0,0,k(),k());
                dmatrix Vxy = V.block(0,k(),k(),m());
                dmatrix Vyx = V.block(k(),0,m(),k());
                dmatrix Vyy = V.block(k(),k(),m(),m());

                dvector o = d.y[i];

                if(Vyy.isZero())
                {
                    res.est[i] = uncertain(Mx, Vxx);
                    res.contrasts[i-1]=0;
                }
                else
                {
                    dmatrix Vinv = pseudoinverse(Vyy); // Vyy.inverse();
                    dvector N = Mx + Vxy * Vinv * (o-My);
                    dmatrix W =Vxx + ((-1) * Vxy) * Vinv * Vyx;
                    assert(isVar(W,"W",i-1,dv(pars)));
                    res.est[i] = uncertain(N,W);

                    if(contrasts && Vyy.determinant() != 0)
                    {
                        vector<double> w = act==ectwls ? d.w[i] : vector<double>();
                        double ll = contrast(act,o,My,Vyy,Vinv, w);
    //                    assert(ll <= 0);
                        res.contrasts[i-1] = ll;
                        res.contrast += ll;
                    }
                }
                g.y.push_back(d.y[i]);
                g.est.push_back(res.est[i].x());
            }
            else
            {
                g.y.push_back(My);
                g.est.push_back(Mx);
            }
            g.z.push_back(d.z[i]);
            res.Ts[i] = this->T(i, pars, g);
            res.Js[i] = this->J(i,pars, g);
            res.Ps[i] = this->P(i,pars, g);
        }
        return res;
    }
private:
    double contrast(const std::vector<double> &v,
                  const data& d, contrasttype act = ectll) const
    {
        evalresult r = eval(v,d,0,ct);
        return r.contrast + contrastaddition(r);
    }
    static double llobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
    {
        seirmodel& f = *((seirmodel*) f_data);
        vector<double> x(f.filter.size());
        unsigned i=0;
        for(unsigned j=0; j<f.filter.size(); j++)
        {
            if(f.filter[j])
                x[j]=f.prs[j];
            else
                x[j]=v[i++];
        }
        double ret = f.contrast(x,f.d,f.ct);
        clog << ret << endl;

        return ret;
    }




    // state variables
    vector<double> prs;
    vector<bool> filter;
    data d;
    contrasttype ct;

public:


    double estimate(const vector<paraminfo>& params,
                    const data& ad,
                    uncertain& res,
                    contrasttype act,
                    vector<bool> omit = vector<bool>(),
                    double timest = 5*60)
    {
        ct = act;
        using namespace nlopt;
        vector<double> upper;
        vector<double> lower;
        vector<double> pars;
        if(omit.size() == 0)
            filter = vector<bool>(params.size(),false);
        else
        {
            assert(omit.size() == params.size());
            filter = omit;
        }
        prs.resize(0);       
        for(unsigned i=0; i<params.size();i++)
        {
            prs.push_back(params[i].initial);
            if(!filter[i])
            {
                upper.push_back(params[i].upper);
                lower.push_back(params[i].lower);
                pars.push_back(params[i].initial);
            }
        }
        d = ad;

        opt o(LN_COBYLA,pars.size());

        o.set_lower_bounds(lower);
        o.set_upper_bounds(upper);

        o.set_ftol_rel(1e-7);
        o.set_maxtime(timest);

        o.set_max_objective(llobj, this);

        nlopt::result oresult;
        double ov;
        try
        {
           oresult = o.optimize(pars, ov);
//oresult = (nlopt::result) 0;
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
        unsigned j=0;
        vector<double> r;
        for(unsigned i=0; i<params.size(); i++)
        {
            if(filter[i])
                r.push_back(params[i].initial);
            else
                r.push_back(pars[j++]);
        }

        try {
            res = uncertain(dv(r), paramvar(r,ad,act));
        } catch (...) {
            res = uncertain(r); // tbd nan
        }

        return ov;
    }



    pair<vector<vector<double>>,vector<vector<double>>>
      sensitivity(const vector<double>& pars,
                        const data& d, contrasttype act ) const
    {
        int numvals = 4;

        dvector dx = 0.001 * dv(pars);
        vector<vector<double>> ret(pars.size(),vector<double>(2 * numvals + 1));
        vector<vector<double>> retx(pars.size(),vector<double>(2 * numvals + 1));
        for(unsigned i=0; i<pars.size(); i++)
        {
            for(int j = -numvals; j <= numvals; j++)
            {
                vector<double> p = pars;
                p[i] += j * dx[i];
                retx[i][j+numvals] = p[i];
                try {
                   evalresult fh = eval(p,d,0,act);
                   ret[i][j+numvals]=fh.contrast;
                } catch (...)
                {
                   ret[i][j+numvals]= std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
        return {retx,ret};
    }





    vector<double> grad(const vector<double>& pars,
                        const data& d, contrasttype act ) const
    {
        evalresult f0 = eval(pars,d,0,act);

        dvector dx =
                sqrt(std::numeric_limits<double>::epsilon()) * dv(pars);
        vector<double> ret(pars.size());
        for(unsigned i=0; i<pars.size(); i++)
        {
            vector<double> p = pars;
            p[i] += dx[i];
            evalresult fh = eval(p,d,0,act);
            ret[i]=(fh.contrast-f0.contrast) / dx[i];
        }
        return ret;
    }

    dmatrix paramvar(const vector<double>& pars,
                    const data& d, contrasttype act) const
    {
        evalresult fixed = eval(pars,d,0,act);
        vector<evalresult> perturbed;
        dvector h = sqrt(std::numeric_limits<double>::epsilon()) * dv(pars);
        for(unsigned i=0; i<h.size(); i++)
            if(h[i] < 1e-10)
                h[i] = 1e-10;
        dvector xph = dv(pars)+h;
        dvector dx = xph - dv(pars);
        for(unsigned i=0; i<pars.size(); i++)
        {
            vector<double> p = pars;
            p[i] = xph[i];
            perturbed.push_back(eval(p,d,0,act));
        }
        dmatrix res(pars.size(),pars.size());
        unsigned t = fixed.contrasts.size();
        for(unsigned i=0; i<t; i++)
        {
            dvector g(perturbed.size());
            for(unsigned j=0; j<perturbed.size(); j++)
                g[j] = (perturbed[j].contrasts[i] - fixed.contrasts[i]) / dx[j];
            res = res + g * g.transpose();            
        }
        try
        {
           return (1.0 / t * res).inverse();
        }
        catch(...)
        {
            cerr << "Cannot invert dmatrix" << endl;
            cerr << 1.0 / t * res;
            throw;
        }
    }

    virtual dmatrix P(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dmatrix J(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual vector<double> I(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual dmatrix Gamma(unsigned /* t */, const vector<double>& /* params */, const G& /*g*/) const
        { dmatrix ret(m(),k()); ret.setZero(); return ret;; }

    virtual double contrastaddition(const evalresult& r) const { return 0; }
    virtual double varfactor(const vector<double>&/* params */) const { return 1; }
    virtual string statelabel(unsigned i) const = 0;
    virtual vector<double> X0(const vector<double>& params, const G& g) const = 0;

    seirmodel(const dmatrix& aH) :  H(aH) {}

    unsigned k() const { return H.cols(); }
    unsigned m() const { return H.rows(); }

protected:
    const dmatrix H;
};

#endif // SEIR_HPP


