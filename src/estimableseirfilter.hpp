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

enum estimationmethod { emwls, emmle, emstdres, emwlsstd, emwlsmle };

template <estimationmethod method=emwls>
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
                res[i]= (g.Y(t,i)-g.Y(former,i))/ (t-former);
        }
        return res;
    }



    virtual double contrast(unsigned t, const vector<double>& pars, G& g, evalparams ep) const
    {
        double mlecoef = method == emwlsmle  ? -ep.additionalcontrastweight : 1.0;
        g.icontrasts[t] = vector<double>(n(),0);
        double s=0;
        if constexpr(method==emwls || method==emwlsstd || method==emwlsmle)
        {
            vector<double> wa = weekave(t,g);
            for(unsigned i=0; i<n(); i++)
            {
                double w = wa[i] < 10 ? 1.0 / 10 : 1.0 / wa[i];
                double x = g.Y(t,i)-
                        (ep.longpredtocontrast
                           ?  g.predlong[t].x()[k()+i]
                          :g.pred[t].x()[k()+i]);
                g.icontrasts[t][i]=w * x * x;
                s += w * x * x;
            }
        }
        if constexpr(method==emmle || method==emwlsmle)
        {
            const uncertain& x = ep.longpredtocontrast
                    ? g.predlong[t]: g.pred[t];
            dmatrix Vyy=x.var().block(k(),k(),n(),n());
            double det = Vyy.determinant();
            if(!det)
            {
                cerr << "Wyy det = 0" << endl;
                cerr << "Params" << endl;
                for(unsigned i=0; i<pars.size(); i++)
                   cerr << pars[i] << endl;
                cerr << m2csv(Vyy) << endl;
                throw( "Vyy det = 0!" );
            }

            dvector d =g.Y(t)-x.x().block(k(),0,n(),1);

            for(unsigned i=0; i<n(); i++)
                g.icontrasts[t][i] += mlecoef * d[i] / sqrt(Vyy(i,i));

            s += mlecoef * (-(n()/2.0) * ::log(2*3.1415926)
                -0.5 * ::log(det)
                -0.5 * (d.transpose()*Vyy.inverse()*d)(0,0));
        }
        return s;
    }

    virtual double contrastaddition(const vector<double>& params, const G& r ,
                                  evalparams ep ) const
    {
        if(method == emstdres || method == emwlsstd)
        {
            unsigned horizon = r.Ysize()-1;
            unsigned esthorizon = horizon - ep.estoffset;

            double contrast = 0;
            for(unsigned i=0; i<n(); i++)
            {
                unsigned ki=k()+i;
                double s=0;
                double s2=0;
                unsigned n=0;
                for(unsigned j=ep.firstcomputedcontrasttime
                     ;j<=esthorizon;j++)
                {
                    const uncertain& x = ep.longpredtocontrast
                            ? r.predlong[j]: r.pred[j];
                    double dd = (r.Y(j)[i]-x.x()[ki]) / sqrt(x.var()(ki,ki));

                    s += dd;
                    s2 += dd * dd;
                    n++;
                }
                double stdev = sqrt( s2 / n - s*s / n / n );
                if(stdev == 0)
                    throw "zero stdev";
                contrast += exp((stdev-1)*(stdev-1));
            }
            return contrast;
        }
        else
            return 0;
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

//        double ret = f.totalcontrast(x,f.d,f.ct);
        clog << r.contrast << endl;
        return r.contrast;
    }

public:


    double estimate(const vector<seirparaminit>& params,
                    const seirdata& ad,
                    const evalparams& ep,
                    uncertain& res,
                    double timest = numeric_limits<double>::infinity()
                    )
    {
        fdata = ad;
        fevalparams = ep;

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

        opt o(LN_COBYLA,nloptpars.size());

        o.set_lower_bounds(lower);
        o.set_upper_bounds(upper);

        o.set_ftol_rel(1e-7);
        o.set_maxtime(timest);

        switch(method)
        {
            case emwls:
            case emwlsstd:
            case emwlsmle:
            case emstdres:
                o.set_min_objective(llobj, this);
            break;
            case emmle:
                o.set_max_objective(llobj, this);
            break;
        default:
            throw "unknown method";
        }

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
