#ifndef SEIR_HPP
#define SEIR_HPP

#include "orpp.hpp"
#include "orpp/matrices.hpp"

using namespace orpp;

inline bool isVar(const matrix& W, const string& label)
{
    for(unsigned j=0; j<W.c();j++)
    {
        if(W(j,j) < -1e-10)
        {
            cerr << label << "("<<j<<","<<j<<")=" << W(j,j) << endl;
            cerr << label << endl;
            cerr << W << endl;
            return false;
        }
    }
    return true;
}

class seirmodel : public object
{
protected:
    struct history { vector<vector<double>> y; vector<vector<double>> z; };
    struct G: public history
    {
        vector<vector<double>> est;
    };

public:
    struct data: public history
    {
        vector<string> dates;
        vector<string> ylabels;
        vector<string> zlabels;
    };
protected:
    matrix T(unsigned t, const vector<double>& pars, const G& g) const
    {
        return P(t,pars,g).transpose()+I(t,pars,g);
    }
    matrix Phi(unsigned t, const vector<double>& pars, unsigned s, const G& g) const
    {
        matrix P = this->P(t,pars,g);
        unsigned k = this->k();
        matrix ret(k,k,0);
        for(unsigned i=0; i<k; i++)
            for(unsigned j=0; j<k; j++)
            {
                ret(i,j) = -P(s,i)*P(s,j);
                if(i==j)
                    ret(i,i)+=P(s,i);
            }
        return ret;
    }

    matrix bigPhi(unsigned t, const vector<double>& pars, unsigned s,const G& g) const
    {
        matrix phi = Phi(t,pars,s,g);
        return block(phi,phi*H.transpose(),H*phi,H*phi*H.transpose());
    }


    matrix bigIwGamma(unsigned t, const vector<double>& pars,
                const vector<double>& N, const G& g ) const
    {
        vector<double> n(N.size());
        for(unsigned i=0; i<N.size(); i++)
            n[i] = max(N[i],0.0);

        matrix ret(k()+m(),k()+m(),0);
        vector<double> Ix = vec(I(t,pars,g) * n);
        matrix D= diag(Ix);
        return block( D, D*H.transpose(), H*D, H*D*H.transpose()
                      + diag(vec(Gamma(t,pars, g)*n)));
    }

    virtual double contrast(const vector<double>& o,const vector<double>& My,
                            const matrix& Vyy, const matrix& Vinv) const
    {
/*        cout << Vyy(0,0) << ", " << vec(o-My)[0] << ", "
             << Vyy(1,1) << ", " << vec(o-My)[1] << ", "
             << ::log(det(Vyy)) << ", " << (transpose(o-My)*Vinv*(o-My))(0,0) << ","
             << 0.5 * ::log(Vyy(0,0)) + 0.5 * vec(o-My)[0] * vec(o-My)[0] / Vyy(0,0) << ","
             << 0.5 * ::log(Vyy(1,1)) + 0.5 * vec(o-My)[1] * vec(o-My)[1] / Vyy(1,1) << ","
             << -(o.size()/2.0) * ::log(2*3.1415926)
                -0.5 * ::log(det(Vyy))
                -0.5 * (transpose(o-My)*Vinv*(o-My))(0,0)
             << endl;*/
//auto oo=o;
//oo[1]=My[1];
//        return -(oo.size()/2.0) * ::log(2*3.1415926)
//            -0.5 * ::log(det(Vyy))
//            -0.5 * (transpose(oo-My)*Vinv*(oo-My))(0,0);
        return -(o.size()/2.0) * ::log(2*3.1415926)
            -0.5 * ::log(det(Vyy))
            -0.5 * (transpose(o-My)*Vinv*(o-My))(0,0);
    }


public:

    struct evalresult
    {
        evalresult(unsigned t, unsigned s, const seirmodel& am) :
            pred(s+1), est(t+1),Ts(s+1),yVs(t+1),
            contrasts(t,0), contrast(0), sm(am)
           { assert(s>=t); }
        vector<uncertain> pred;
        vector<uncertain> est;
        vector<matrix> Ts;
        vector<matrix> yVs;
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
                    for(unsigned i=0; i<d.y[s].size(); i++)
                            o<< d.y[s][i] << ",";
                    for(unsigned i=0; i<d.z[s].size(); i++)
                            o<< d.z[s][i] << ",";
                }
                o << endl;
            }
        }
    };

/*    pair<vector<vector<double>>,vector<vector<double>>> evaltrend(const vector<double>& pars, unsigned ds=0) const
    {
        unsigned t = obs().size();
        vector<vector<double>> resobs(t+ds);
        vector<vector<double>> resx(t+ds);
        vector<double> x0 = X0(pars);
        vector<double> lastx = x0;
        for(unsigned i=0; i<t+ds; i++)
        {
            matrix T = this->T(i, pars);
            vector<double> x = vec(T*lastx);

            matrix THT = stackv(T,H*T);
            resx[i] = x;
            resobs[i] = vec(H*x);
            lastx = x;
        }
        return {resx, resobs};
    }*/

    template <bool contrasts=true>
    evalresult eval(const vector<double>& pars,
                    const data& d, unsigned ds=0) const
    {
        assert(d.y.size());
        unsigned t = d.y.size()-1;
        assert(ds >= 0);
        assert(d.z.size()>=t+ds);
        evalresult res(t,t+ds,*this);
        G g;
        g.y = {d.y[0]};
        g.z = {d.z[0]};
        vector<double> x0 = X0(pars,g);
        res.pred[0]=stackv(x0, vec(H*x0));
        g.est = {x0};
        res.est[0]=x0;
        res.Ts[0]=this->T(0, pars, g);
        unsigned i=1;
        for(; i<= t+ds; i++)
        {
            matrix& T = res.Ts[i-1];
            matrix THT = stackv(T,H*T);
            matrix UhU = stackv(U,H*U);

            vector<double> N;
            matrix W;
            if(i<=t+1)
            {
                N = res.est[i-1].x();
                W = res.est[i-1].var();
                assert(isVar(W,"W1"));
            }
            else
            {
                N = res.pred[i-1].x();
                N.resize(k());
                W = res.pred[i-1].var().submatrix(0,0,k(),k());
                assert(isVar(W,"W2"));
            }
            vector<double> Mx = vec(T * N + U * g.z[i-1]);
            vector<double> My = vec(H * Mx);

            matrix V = THT * W * THT.transpose();
            assert(isVar(V,"V2"));
            for(unsigned j=0; j<k(); j++)
                V = V + max(N[j],0.0) * bigPhi(i-1,pars,j,g);
            assert(isVar(V,"V3"));

            V = V + bigIwGamma(i-1, pars, N,g);
            assert(isVar(V,"V4"));

            res.pred[i] = uncertain(stackv(Mx, My),V);
            if(i<=t)
            {
                matrix Vxx = V.submatrix(0,0,k(),k());
                matrix Vxy = V.submatrix(0,k(),k(),m());
                matrix Vyx = V.submatrix(k(),0,m(),k());
                matrix Vyy = V.submatrix(k(),k(),m(),m());

                vector<double> o = d.y[i];

                if(Vyy.iszero())
                {
                    res.est[i] = uncertain(Mx, Vxx);
                    res.contrasts[i-1]=0;
                }
                else
                {
                    matrix Vinv = pseudoinverse(Vyy); // Vyy.inverse();
                    vector<double> N = vec(Mx + Vxy * Vinv * (o-My));
                    matrix W =Vxx + ((-1) * Vxy) * Vinv * Vyx;
                    assert(isVar(W,"W"));
                    res.est[i] = uncertain(N,W);

                    if(contrasts && det(Vyy) != 0)
                    {
                        double ll = contrast(o,My,Vyy,Vinv);
                        res.yVs[i] = Vinv;
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
        }
        return res;
    }
private:
    double contrast(const std::vector<double> &v,
                  const data& d) const
    {
        evalresult r = eval(v,d);        
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
        double ret = f.contrast(x,f.d);
        clog << ret << endl;

        return ret;
    }

/*    static double trendobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
    {
        seirmodel& f = *((seirmodel*) f_data);
        vector<double> x(f.filter.size());
        unsigned i=0;
        for(unsigned j=0; j<f.filter.size(); j++)
            if(f.filter[j])
                x[j]=f.prs[j];
            else
                x[j]=v[i++];
        vector<vector<double>> r = f.evaltrend(x).second;
        vector<vector<double>> os = f.obs();
        double ret = 0;
        for(unsigned i=0; i<os.size(); i++)
        {
            for(unsigned j=0; j<os[i].size(); j++)
                ret += pow(os[i][j]-r[i][j],2);
        }
        clog << ret << endl;
        return ret;
    }*/



    // state variables
    vector<double> prs;
    vector<bool> filter;
    data d;
public:


    double estimate(const vector<paraminfo>& params,
                    const data& ad,
                    uncertain& res,
                    vector<bool> omit = vector<bool>())
    {
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
        o.set_maxtime(6*60);

        o.set_max_objective(llobj, this);

        nlopt::result oresult;
        double ov;
        try
        {
            oresult = o.optimize(pars, ov);
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
            res = uncertain(r, paramvar(r,ad));
        } catch (...) {
            res = uncertain(r); // tbd nan
        }

        return ov;
    }


    template <bool negative = false>
    vector<double> grad(const vector<double>& pars,
                        const data& d) const
    {
        evalresult f0 = eval(pars,d);

        vector<double> h = sqrt(std::numeric_limits<double>::epsilon()) * pars;
        vector<double> xph = negative ? pars-h : pars+h;
        vector<double> dx = negative ? pars-xph : xph - pars;
        vector<double> ret(pars.size());
        for(unsigned i=0; i<pars.size(); i++)
        {
            vector<double> p = pars;
            p[i] = xph[i];
            evalresult fh = eval(p,d);
            ret[i]=(fh.contrast-f0.contrast) / dx[i];
        }
        return ret;
    }

    matrix paramvar(const vector<double>& pars,
                    const data& d) const
    {
        evalresult fixed = eval(pars,d);
        vector<evalresult> perturbed;
        vector<double> h = sqrt(std::numeric_limits<double>::epsilon()) * pars;
        for(unsigned i=0; i<h.size(); i++)
            if(h[i] < 1e-10)
                h[i] = 1e-10;
        vector<double> xph = pars+h;
        vector<double> dx = xph - pars;
        for(unsigned i=0; i<pars.size(); i++)
        {
            vector<double> p = pars;
            p[i] = xph[i];
            perturbed.push_back(eval(p,d));
        }
        matrix res(pars.size(),pars.size());
        unsigned t = fixed.contrasts.size();
        for(unsigned i=0; i<t; i++)
        {
            vector<double> g(perturbed.size());
            for(unsigned j=0; j<perturbed.size(); j++)
                g[j] = (perturbed[j].contrasts[i] - fixed.contrasts[i]) / dx[j];
            res = res + g * transpose(g);
        }
        try
        {
           return (1.0 / t * res).inverse();
        }
        catch(...)
        {
            cerr << "Cannot invert matrix" << endl;
            cerr << 1.0 / t * res;
            throw;
        }
    }

    virtual matrix P(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual matrix I(unsigned t, const vector<double>& params, const G& g) const = 0;
    virtual matrix Gamma(unsigned /* t */, const vector<double>& /* params */, const G& /*g*/) const
        { return matrix(m(),k(),0); }

    virtual double contrastaddition(const evalresult& r) const { return 0; }
    virtual string statelabel(unsigned i) const = 0;
    virtual vector<double> X0(const vector<double>& params, const G& g) const = 0;

    seirmodel(const matrix& aH, const matrix& aU) :  H(aH), U(aU) {}

    unsigned k() const { return H.c(); }
    unsigned m() const { return H.r(); }

protected:
    const matrix H;
    const matrix U;
};

#endif // SEIR_HPP


