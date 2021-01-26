#ifndef ESTIMABLESEIRFILTER_HPP
#define ESTIMABLESEIRFILTER_HPP

#include "orpp/estimation.hpp"
#include "seirfilter.hpp"

using namespace orpp;
using namespace std;



struct seirparaminit
{
    std::string name;
    double initial;
    double lower;
    double upper;
    bool omit;
};

template <bool wls=true>
class estimableseirfilter : virtual public seirfilter
{

public:

    estimableseirfilter()
    {
    }

    vector<double> weekave(unsigned t, const G& g) const
    {
        vector<double> res(g.n(),0);
        if(t>0)
        {
            unsigned former = t>7 ? t-7 : 0;
            for(unsigned i=0; i<n(); i++)
                res[i]= (g.y[t][i]-g.y[former][i])/ (t-former);
        }
        return res;
    }

    virtual double contrast(unsigned t, const vector<double>&, G& g, bool longpredtocontrast) const
    {
        if constexpr(wls)
        {
            vector<double> wa = weekave(t,g);
            double s=0;
            vector<double> ics(n(),0);
            for(unsigned i=0; i<n(); i++)
            {
                double w = wa[i] < 10 ? 1.0 / 10 : 1.0 / wa[i];
                double x = g.y[t][i]-
                        (longpredtocontrast
                           ?  g.predlong[t].x()[k()+i]
                          :g.pred[t].x()[k()+i]);
                ics[i]=w * x * x;
                s += w * x * x;
            }
            g.icontrasts[t] = ics;
            return s;
        }
        else
        {
            const uncertain& x = longpredtocontrast
                    ? g.predlong[t]: g.pred[t];
            dmatrix Vyy=x.var().block(k(),k(),n(),n());
            double det = Vyy.determinant();
            if(!det)
            {
                cerr << m2csv(Vyy) << endl;
                throw( "Vyy det = 0!" );
            }

            dvector d =g.y[t]-x.x().block(k(),0,n(),1);
            vector<double> ics(n(),0);

            for(unsigned i=0; i<n(); i++)
            {
                ics[i] = d[i] / sqrt(Vyy(i,i));
            }

            g.icontrasts[t] = ics;

            return -(n()/2.0) * ::log(2*3.1415926)
                -0.5 * ::log(det)
                -0.5 * (d.transpose()*Vyy.inverse()*d)(0,0);
        }
    }
private:

    // state variables
    vector<double> fparams;
    vector<bool> ffilter;
    evalparams fevalparams;
    seirdata fdata;


    static double llobj(const std::vector<double> &v, std::vector<double> &, void* f_data)
    {
        estimableseirfilter& f = *((estimableseirfilter*) f_data);
        vector<double> x(f.ffilter.size());
        unsigned i=0;
        for(unsigned j=0; j<f.ffilter.size(); j++)
            {
            if(f.ffilter[j])
                x[j]=f.fparams[j];
            else
                x[j]=v[i++];
        }

        G r = f.eval(x,f.fdata, f.fevalparams);
        double ret = r.contrast+f.contrastaddition(v,r);

//        double ret = f.totalcontrast(x,f.d,f.ct);
        clog << ret << endl;
        return ret;
    }

public:


    double estimate(const vector<seirparaminit>& params,
                    const seirdata& ad,
                    const evalparams& ep,
                    uncertain& res,
                    double timest = numeric_limits<double>::infinity()
                    )
    {
        using namespace nlopt;
        vector<double> upper;
        vector<double> lower;
        vector<double> nloptpars;
        fparams.resize(0);
        ffilter.resize(params.size());
        for(unsigned i=0; i<params.size();i++)
        {
            fparams.push_back(params[i].initial);
            ffilter[i] = params[i].omit;
            if(!ffilter[i])
            {
                upper.push_back(params[i].upper);
                lower.push_back(params[i].lower);
                nloptpars.push_back(params[i].initial);
            }
        }
        fdata = ad;
        fevalparams = ep;

        opt o(LN_COBYLA,nloptpars.size());

        o.set_lower_bounds(lower);
        o.set_upper_bounds(upper);

        o.set_ftol_rel(1e-7);
        o.set_maxtime(timest);

        if(wls)
            o.set_min_objective(llobj, this);
        else
            o.set_max_objective(llobj, this);

        nlopt::result oresult;
        double ov;
        try
        {
           oresult = o.optimize(nloptpars, ov);
        }
        catch(std::exception &e)
        {
            throw std::runtime_error(e.what());
        }
        if(oresult < 0)
        {
            ostringstream es;
            es << "Negative result of nlopt: " << oresult << endl;
            throw std::runtime_error(es.str());
        }
        unsigned j=0;
        vector<double> r;
        for(unsigned i=0; i<params.size(); i++)
        {
            if(ffilter[i])
                r.push_back(params[i].initial);
            else
                r.push_back(nloptpars[j++]);
        }

        G fixed = eval(r,ad, ep);
        vector<G> perturbed;
        dvector h = sqrt(std::numeric_limits<double>::epsilon()) * dv(r);
        for(unsigned i=0; i<h.size(); i++)
            if(h[i] < 1e-10)
                h[i] = 1e-10;
        dvector xph = dv(r)+h;
        dvector dx = xph - dv(r);
        for(unsigned i=0; i<r.size(); i++)
        {
            vector<double> p = r;
            p[i] = xph[i];
            perturbed.push_back(eval(p,ad,ep));
        }
        dmatrix x(r.size(),r.size());
        x.setZero();
        unsigned t = fixed.contrasts.size();
        for(unsigned i=0; i<t; i++)
        {
            dvector g(perturbed.size());
            for(unsigned j=0; j<perturbed.size(); j++)
                g[j] = (perturbed[j].contrasts[i] - fixed.contrasts[i]) / dx[j];
            x += g * g.transpose();
        }
        try
        {
           res = uncertain(dv(r),(1.0 / t * x).inverse());
        }
        catch(...)
        {
            cerr << "Cannot invert dmatrix" << endl;
            cerr << m2csv(1.0 / t * x);
            res = uncertain(r);
        }
        return ov;
    }
};


#endif // ESTIMABLESEIRFILTER_HPP
