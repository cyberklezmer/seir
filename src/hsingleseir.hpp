#ifndef HSINGLESEIR_HPP
#define HSINGLESEIR_HPP

#include "hseir.hpp"

class hsingleseir: public hcohortseir
{
public:
    enum obscolums { /* R */ RA, RS, H, G, DH, DO,  numobscolumns };

    hsingleseir() : hcohortseir(1,{ 10699000 }) {}

    virtual unsigned n() const
    {
        return numobscolumns;
    }
    virtual dmatrix F(unsigned, const vector<double>& , const struct G& ) const
    {

//     E,  Ia,  Ip,  Is,  Iu  R, Edelta, Iadelta, Ipdelta,  Isdelta,  Rdelta,  Hdelta, Rhdelta, Ddelta, Dhdelta,  Isd, Rd, Hd, Rhd, Dd, Dhd
       dmatrix ret(n(),k());
       ret <<
/*RA*/   0,  0,   0,   0,   0,  0, 0,      1,       1,        1,        1,       1,      1,       1,      1,        0,   0,  0,  0,   0,  0,
/*RS*/     0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      0,       0,      0,        1,   1,  1,  1,   1,  1,
/*R  0,  0,   0,   0,   0,  0, 0,      1,       1,        1,        1,       1,      1,       1,      1,        1,   1,  1,  1,   1,  1,*/
/*H*/  0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       1,      1,       0,      1,        0,   0,  1,  1,   0,  1,
/*G*/  0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      1,       0,      1,        0,   0,  0,  1,   0,  1,
/*DH*/ 0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      0,       0,      1,        0,   0,  0,  0,   0,  1,
/*DO*/ 0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      0,       1,      0,        0,   0,  0,  0,   1,  0;
        return ret;
    }
    virtual dvector I(unsigned t, const vector<double>&, const struct G& g) const
    {
        dvector ret(k());
        ret.setZero();
        double s=0;
        for(unsigned i=I0; i<=I80; i++)
            s += g.Z(t,i);
        ret[hpartial::E] = s;
        return ret;
    }

    virtual dvector gamma(unsigned /* t */, const vector<double>& params, const struct G& /*g*/) const
    {
        dvector ret(n());
        ret.setZero();
        ret[RA] = params[rcoef];
        ret[RS] = params[rcoef];
        //ret[RS] = params[rscoef];
        ret[G] = params[gcoef];
        ret[H] = params[gcoef];
        ret[DH] = ret[DO] = params[dcoef];
        return ret;
    }

};

class singledatareader : public datareader
{
    virtual vector<string> obslabels() const
    {
        return   { "RA", "RS",  "H", "G", "DH", "DO"};
    }

    virtual unsigned numobs() const
    {
        return hsingleseir::numobscolumns;
    }
    virtual vector<double> obs2obs(const vector<double>& src) const
    {
        vector<double> dst(hsingleseir::numobscolumns);

        dst[hsingleseir::RA] = src[RA];
        dst[hsingleseir::RS] = src[RS];
//        dst[hsingleseir::R] = src[RS]+src[RA];
        dst[hsingleseir::H] = src[HY] + src[HO];
        dst[hsingleseir::G] = src[DHY] + src[DHO]
                                      + src[RHY] + src[RHO];
        dst[hsingleseir::DH] = src[DHY] + src[DHO];
        dst[hsingleseir::DO] = src[DOY] + src[DOO];
        return dst;
    }
};


#endif // HSINGLESEIR_HPP
