#ifndef COHORTSEIR_HPP
#define COHORTSEIR_HPP


#include "seirfilter.hpp"

class partialseir
{
public:
   virtual dmatrix P(const vector<double>& params) const = 0;
   virtual dmatrix hatB(const vector<double>& params) const = 0;
   virtual unsigned k() const = 0;

   virtual std::string statelabel(unsigned i) const
   {
       return std::to_string(i);
   }

   virtual double susceptible( double /* population */, const dvector& /* compartments */) const
   {
       return 1;
   }
};

template <typename C>
class cohortseir : virtual public seirfilter
{
public:
     cohortseir(unsigned numcohorts, const vector<double>& population) :
         fnumcohorts(numcohorts), fpopulation(population)
     {
         assert(fpopulation.size()==fnumcohorts);
     }

     virtual dmatrix mixingmatrix(unsigned, const vector<double>& , const G& ) const
     {
         dmatrix ret(fnumcohorts,fnumcohorts);
         ret.setIdentity();
         return ret;
     }

     virtual std::vector<std::vector<double>> pars2pars(
             unsigned t, const vector<double>& params, const G& g) const = 0;

     virtual std::string statelabel(unsigned i) const
     {
         assert(i < numcohorts() * fpartial.k());
         unsigned cn = i / fpartial.k();
         unsigned sn = i % fpartial.k();
         return cohortlabel(cn) + "_" + fpartial.statelabel(sn);
     }

     virtual std::string cohortlabel(unsigned i) const
     {
         assert(i<numcohorts());
         return std::to_string(i);
     }
     unsigned numcohorts() const { return fnumcohorts; }

     virtual dmatrix P(unsigned t, const vector<double>& params, const G& g) const
     {
         dmatrix ret(k(),k());
         ret.setZero();
         unsigned pk = fpartial.k();
         auto pars = pars2pars(t,params,g);
         for(unsigned i=0; i<numcohorts() ; i++)
             ret.block(i*pk,i*pk,pk,pk) = fpartial.P(pars[i]);
         return ret;
     }
     virtual dmatrix hatB(unsigned t, const vector<double>& params, const G& g) const
     {
         dmatrix ret(k(),k());
         ret.setZero();
         unsigned nc = numcohorts();
         unsigned pk = fpartial.k();
         auto pars = pars2pars(t,params,g);
         dmatrix mm = mixingmatrix(t,params,g);
         const dvector& c = g.est[t].x();
         for(unsigned i=0; i<nc; i++)
         {
             auto b = fpartial.hatB(pars[i]);
             double s = fpartial.susceptible(fpopulation[i],
                           c.block(i*pk,0,pk,1));
             for(unsigned j=0; j<nc; j++)
                 ret.block(i*pk,j*pk,pk,pk)=mm(i,j)*b*s;
         }
         return ret;
     }

     virtual unsigned k() const
     {
         return numcohorts()*fpartial.k();
     }

     const C& partial() const { return fpartial; }

     virtual double abssusceptible(unsigned t, const vector<double>& params, const G& g) const
     {
         unsigned nc = numcohorts();
         unsigned pk = fpartial.k();
         auto pars = pars2pars(t,params,g);
         const dvector& c = g.est[t].x();
         double s = 0;
         for(unsigned i=0; i<nc; i++)
             s += fpopulation[i] *  fpartial.susceptible(fpopulation[i],
                           c.block(i*pk,0,pk,1));
         return s;
     }

     virtual dmatrix activesubmatrix(const dmatrix &T) const
     {
         vector<vector<dmatrix>> s(fnumcohorts,vector<dmatrix>(fnumcohorts));
         unsigned pk = fpartial.k();
         for(unsigned i=0; i<fnumcohorts; i++)
             for(unsigned j=0; j<fnumcohorts; j++)
                 s[i][j] = fpartial.activesubmatrix(T.block(i*pk,j*pk,pk,pk));
         unsigned pm = s[0][0].cols();
         dmatrix ret(fnumcohorts*pm,fnumcohorts*pm);
         ret.setZero();

         for(unsigned i=0; i<fnumcohorts; i++)
             for(unsigned j=0; j<fnumcohorts; j++)
             {
                 assert(s[i][j].cols() == pm);
                 assert(s[i][j].rows() == pm);
                 ret.block(pm * i, pm * j, pm,pm) = s[i][j];
             }
         return ret;
     }
private:
     unsigned fnumcohorts;
     vector<double> fpopulation;;
     C fpartial;
};


#endif // COHORTSEIR_HPP
