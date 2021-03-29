#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <vector>
#include <cmath>
#include <assert.h>

using namespace std;

class gauss
{
public:
    const double scale;
    const double mean;
    gauss(double amean=0, double ascale = 1.0) : scale(ascale), mean(amean) { }
    double operator () (const double& x) const
    {
       return std::erfc(-(x-mean)/scale/std::sqrt(2))/2;
    }
};

class mixfunction
{

    static vector<gauss> makefns(const vector<double>& means, const vector<double>& scales)
    {
        vector<gauss> res;
        for(unsigned i=1; i<means.size(); i++)
            res.push_back(gauss(means[i],scales[i]));
        return res;
    }

    vector<gauss> basefns;

public:
    mixfunction(const vector<double>& means, const vector<double>& scales):
        basefns(makefns(means,scales))
    {}

    double operator () (double x) const
    {
        double s=0;
        double p=1.0 /basefns.size();

        for(unsigned i=0; i<basefns.size(); i++)
            s+=p * basefns[i](x);
        return s;
    }

    double inv(double y) const
    {
        assert(y>0);
        assert(y<1);
        double lo = -1.0;
        double hi = 1.0;
        for(;;)
        {
            if((*this)(hi) > y)
                break;
            hi *= 2.0;
        }
        for(;;)
        {
            if((*this)(lo) < y)
                break;
            lo *= 2.0;
        }
        for(;;)
        {
            double mid = (hi+lo) / 2.0;
            double fx = (*this)(mid);
            if(fx < y)
                lo = mid;
            else
                hi = mid;
            if(hi-lo < 0.00001)
                return mid;
        }
        assert(0);
        return 0.0;
    }
};


#endif // FUNCTIONS_HPP
