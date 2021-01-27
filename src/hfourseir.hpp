#ifndef HFOURSEIR_HPP
#define HFOURSEIR_HPP

#include "hseir.hpp"

class hfourseir: public hcohortseir
{

public:
    enum obscolums { /* R */
       RA, RS,
       HY, HO, GY, GO, DHY, DHO, // DOY, DOO, can  computed from others
       D0 , D20, D65, D80, R0, R20, R65, /* R80 linearly dependent */
       numobscolumns };

    static dmatrix makeF()
    {
       dmatrix ret(numobscolumns ,hpartial::numstates * 4);
       ret.setZero();

       for(unsigned i=0; i<4; i++)
       {
           auto pk = hpartial::numstates;
           unsigned offset = i*pk;

           using row = Eigen::Matrix<double,1,hpartial::numstates>;

           //     E,  Ia,  Ip,  Is,  Iu,  R, Edelta, Iadelta, Ipdelta,  Isdelta,  Rdelta,  Hdelta, Rhdelta, Ddelta, Dhdelta,  Isd, Rd, Hd, Rhd, Dd, Dhd


           ret.block(RA,offset,1,pk) = row(
                {  0,  0,   0,   0,   0,  0, 0,      1,       1,        1,        1,       1,      1,       1,      1,        0,   0,  0,  0,   0,  0}
             );

           ret.block(RS,offset,1,pk) = row(
                {  0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      0,       0,      0,        1,   1,  1,  1,   1,  1}
             );

           unsigned oi = i >= 2;

           ret.block(HY+oi,offset,1,pk) = row(
           /*H*/ {  0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       1,      1,       0,      1,        0,   0,  1,  1,   0,  1}
             );

           ret.block(GY+oi,offset,1,pk) = row(
           /*G*/ {  0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      1,       0,      1,        0,   0,  0,  1,   0,  1}
             );

           ret.block(DHY+oi,offset,1,pk) = row(
           /*DH*/ { 0,  0,   0,   0,   0,  0, 0,      0,       0,        0,        0,       0,      0,       0,      1,        0,   0,  0,  0,   0,  1}
             );

           //     E,  Ia,  Ip,  Is,  Iu,  R, Edelta, Iadelta, Ipdelta,  Isdelta,  Rdelta,  Hdelta, Rhdelta, Ddelta, Dhdelta,  Isd, Rd, Hd, Rhd, Dd, Dhd

           ret.block(D0+i,offset,1,pk) = row(
           /*D*/ {0,  0,   0,   0,   0,   0,  0,     0,       0,        0,        0,       0,      0,       1,      1,        0,   0,  0,  0,   1,  1}
             );
           if(i<3)
               ret.block(R0+i,offset,1,pk) = row(
               /*R*/ {0,  0,   0,   0,   0,   0,  0,     1,       1,        1,        1,       1,      1,       1,      1,        1,   1,  1,  1,   1,  1}
                 );
       }

       Eigen::FullPivLU<dmatrix> lu_decomp(ret);
       auto rank = lu_decomp.rank();
       if(rank!=numobscolumns)
       {
           cout << m2csv(ret) << endl;
           throw "bad rank of F";
       }
       return ret;
    }


    hfourseir() : hcohortseir(4,{2188232,6374077,1690530,441100}),
        fF(makeF())
       {}

    virtual unsigned n() const
    {
        return numobscolumns;
    }

    virtual dmatrix mixingmatrix(unsigned, const vector<double>& , const G& ) const
    {
        dmatrix ret(4,4);
        ret <<
               0.544423361025132,	0.108333151559267,	0.0212598890307393,	0.0212598890307393,
               0.602958310237031,	0.918820576671648,	0.119606940083756,	0.119606940083756,
               0.08979754354593,	0.0975146549418154,	0.140321057525071,	0.140321057525071,
               0.0240655902239016,	0.0179033162039532,	0.0153248487604178,	0.0587329779527297;

         return ret;
    }


    virtual dmatrix F(unsigned, const vector<double>& , const struct G& ) const
    {
       return fF;
    }

    virtual dvector I(unsigned t, const vector<double>&, const struct G& g) const
    {
        dvector ret(k());
        ret.setZero();
        unsigned j=0;
        for(unsigned i=I0; i<=I80; i++,j++)
            ret[hpartial::E + j*partial().k()] = g.z[t][i];
        return ret;
    }

    virtual double hats(const vector<double>&/* params */) const
    {
        return 10;
    }

    virtual dmatrix Gamma(unsigned /* t */, const vector<double>& /* params */, const struct G& ) const
    {
        dmatrix ret(n(),k());
        ret.setZero();
return ret;
        for(unsigned i=0; i<4; i++)
        {
            unsigned offset = i * partial().k();
            ret(RA,offset+hpartial::E)=ret(RA,offset+hpartial::Ia)
                =ret(RA,offset + hpartial::Ia) = 0.16;
            ret(RS,offset+hpartial::Is)=ret(RS,offset+hpartial::Iu)=1.45;

            unsigned oi = i >= 2;

            ret(HY+oi,offset+hpartial::Hd)=ret(HY+oi,offset+hpartial::Hdelta)=0.19;
            ret(DHY+oi,offset+hpartial::Hd)=ret(DHY+oi,offset+hpartial::Hdelta)=0.0008;
            if(i<3)
                ret(R0+i,offset +hpartial::E)=ret(R0+i,offset+hpartial::Ia)
                    =ret(R0+i,offset + hpartial::Ia) = 0.16;

        }
        return ret;
    }
private:
    dmatrix fF;
};


class fourdatareader : public datareader
{
    virtual vector<string> obslabels() const
    {
        return   { "RA", "RS",
            "HY", "HO","GY","GO","DHY","DHO",
            "D0","D20","D65","D80", "R0","R20","R65" };
    }

    virtual unsigned numobs() const
    {
        return hfourseir::numobscolumns;
    }
    virtual vector<double> obs2obs(const vector<double>& src) const
    {
        vector<double> dst(numobs());

         dst[hfourseir::RA] = src[RA];
        dst[hfourseir::RS] = src[RS];
        dst[hfourseir::HY] = src[HY];
        dst[hfourseir::HO] = src[HO];
        dst[hfourseir::GY] = src[DHY] + src[RHY];
        dst[hfourseir::GO] = src[DHO] + src[RHO];
        dst[hfourseir::DHY] = src[DHY];
        dst[hfourseir::DHO] = src[DHO];
        dst[hfourseir::D0] = src[D0];
        dst[hfourseir::D20] = src[D20];
        dst[hfourseir::D65] = src[D65];
        dst[hfourseir::D80] = src[D80];
        dst[hfourseir::R0] = src[R0];
        dst[hfourseir::R20] = src[R20];
        dst[hfourseir::R65] = src[R65];

        return dst;
    }
};

#endif // HFOURSEIR_HPP
