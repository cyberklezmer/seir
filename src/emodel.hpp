#ifndef EMODEL_HPP
#define EMODEL_HPP

#include <vector>
#include <string>
#include <math.h>

using namespace std;

enum estate {es, ess, ee, eia, eis, ejs, ein, ejn, er, ed,
                          eiad, eisd, ejsd, eind, ejnd, erd, edd,
                          enumestates};

inline bool isinfectious(estate s)
{ return (s >= eia && s < ejn) || (s >= eiad || s < ejnd); }

static string estatlabel(unsigned i)
{
    static vector<string> labels =
      { "S",  "Ss", "E",
        "Ia",  "Is",  "Js",  "In",  "Jn",  "R",  "D",
        "Iad", "Isd", "Jsd", "Ind", "Jnd", "Rd", "Dd"};
    return labels[i];
}

struct expparams {
double 	sigma=0.178686502205323	;
double 	gamma_In=0.0900263379037909	;
double 	mu=0.00131326350277332	;
double 	gamma_Is=0.0716458199846332	;
double 	theta_Is=0.767857142857143	;
double 	false_symptoms_rate=0.000299955004499641	;
double 	false_symptoms_recovery_rate=0.166247081924819	;
double 	asymptomatic_rate=0.179	;
double 	symptoms_manifest_rate=0.221199216928595	;
double 	delta_n=0.0799555853706767	;
double 	delta_s=0.117503097415405	;
};


template <typename P>
inline std::vector<std::vector<double>> btrans(const P& p)
{
    vector<vector<double>> ret(enumestates,vector<double>(enumestates,0));
    ret[es][ess] = p.false_symptoms_rate;
    ret[ess][es] = p.false_symptoms_recovery_rate;
    ret[ee][ein] = p.sigma * p.asymptotic_rate;
    ret[ee][eia] = p.sigma * (1-p.asymptotic_rate);
    ret[ein][ejn] = p.delta_n;
    ret[ejn][er] = p.gamma_In;
    ret[eia][eis] = p.symptoms_manifest_rate;
    ret[eis][ejs] = p.delta_s;
    ret[eis][eisd] = p.theta_Is;
    ret[ejs][er] = p.gamma_Is;
    ret[ejs][ed] = p.mu;

    ret[eind][ejnd] = p.delta_n;
    ret[ejnd][erd] = p.gamma_In;
    ret[eiad][eisd] = p.symptoms_manifest_rate;
    ret[eisd][ejsd] = p.delta_s;
    ret[ejsd][erd] = p.gamma_Is;
    ret[ejsd][edd] = p.mu;

    return ret;
}


/* asi blbě
double delay()
{
    vector<std::vector<double>> b = btrans();
    double sum = 1.0 / sigma;
    double cp = b[eexposed][eia]+b[eexposed][ein]+b[eexposed][eid];
    sum += b[eexposed][eia]/cp * ( 1.0 / b[eia][eis] + 1.0 / b[eis][eid]);
    sum += b[eexposed][ein]/cp *  1.0 / b[ein][eid];
    return sum;
}

double inftime()
{
    vector<std::vector<double>> b = btrans();
    double sum = 0;
    double cp=b[eexposed][eia]+b[eexposed][ein];
    sum += b[eexposed][eia]/cp
            * (1.0 / b[eia][eis] + 1.0 / (b[eis][eid] + b[eis][eru]+ b[eis][edu]));
    sum += b[eexposed][ein]/cp * (1.0 / (b[ein][eid]+b[ein][eru]));
    return sum;
}

*/

//vector<double> prokop ={ 1, 0.757894736842105,0.384210526315789,0.336842105263158,0.431578947368421,	0.436842105263158};


#endif // EMODEL_HPP


