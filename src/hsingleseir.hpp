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
            s += g.z[t][i];
        ret[hpartial::E] = s;
        return ret;
    }

    virtual double hats(const vector<double>&/* params */) const
    {
        return 10;
    }


    virtual dmatrix Gamma(unsigned /* t */, const vector<double>&  params , const struct G& ) const
    {
        dmatrix ret(n(),k());
        ret.setZero();
return ret;
        ret(RA,hpartial::E)=ret(RA,hpartial::Ia)=ret(RA,hpartial::Ia)=0.16;
        ret(RS,hpartial::Is)=ret(RS,hpartial::Iu)=1.45;
        ret(H,hpartial::Hd)=ret(H,hpartial::Hdelta)=0.19;
        ret(DH,hpartial::Hd)=ret(DH,hpartial::Hdelta)=0.0008;
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
