#ifndef CCTOOLS_HPP
#define CCTOOLS_HPP

#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <assert.h>
#include <math.h>
#include <time.h>


using namespace std;
namespace cctools
{

class sys
{
    static sys& self()
    {
       static sys s;
       return s;
    }
public:
    sys() : fout(0), flog(0), ferr(0), funiform(0.0,1.0)
    {
    }
    ~sys()
    {
    }

    /// Sets an existing stream as the text output of the library.
    static void setout(std::ostream& o)
    {
        self().fout = &o;
    }

    /// Resets the standard ouptut of the library,
    /// redirecting it back to \p std::cout
    static void resetout()
    {
        self().fout = 0;
    }

    static void setlog(std::ostream& l)
    {
        self().flog = &l;
    }
    static void resetlog()
    {
        self().flog=0;
    }

    static void seterr(std::ostream& e)
    {
        self().ferr = &e;
    }
    static void reseterr()
    {
        self().ferr = 0;
    }

    static std::ostream& out()
    {
        return self().fout ? *(self().fout) : std::cout;
    }
    static std::ostream& log()
    {
        return self().flog ? *(self().flog) : std::clog;
    }
    static std::ostream& err()
    {
        return self().ferr ? *(self().ferr) : std::cerr;
    }

    /// Resets the seed of the random generator according to computer time.
    static void seed()
    {
        self().fengine.seed(time(0));
    }

    /// Resets the seed of the random generator.
    static void seed(unsigned int aseed)
    {
        self().fengine.seed(aseed);
    }

    /// Returns a pesudo-random observaton from the uniform distribution on [0,1]
    static double uniform()
    {
        return self().funiform(self().fengine);
    }

private:
    std::ostream* fout;
    std::ostream* flog;
    std::ostream* ferr;
    std::default_random_engine fengine;
    std::uniform_real_distribution<double> funiform;
};

inline double pi() { return std::atan(1)*4; }

template <char sep>
inline void csvout(ostream& os, const string& s)
{
    os << '"';
    for(unsigned i=0; i<s.length(); i++)
        if(s[i] == '"')
            os << '"' << '"';
        else
            os << s[i];
    os << '"' << sep;
}

inline unsigned countfalses(const vector<bool>& m)
{
    unsigned n=0;
    for(unsigned i=0; i<m.size(); i++)
        if(!m[i])
            n++;
    return n;
}

template <char sep>
class csv : private vector<vector<string>>
{
    static string getstr(istream& is)
    {
        char c;
        is >> c;
        bool quotes = false;
        if(c=='"')
        {
            quotes = true;
            is >> c;
        }

        string r;
        for(;;)
        {
            if((signed)c > 128)
                throw "Czech characters cannot be used";
            if((c==sep && !quotes) || is.eof())
                break;
            if(c=='"' && quotes)
            {
               is >> c;
               if(c!='"')
                  break;
            }
            r.push_back(c);
            is >> c;
        }
        return r;
    }

    static vector<vector<string>> getvecs(const string& fn)
    {
        vector<vector<string>> res;
        ifstream file(fn);
        if(!file.is_open())
        {
            ostringstream es;
            es << "cannot open csv file: " << fn ;
            throw es.str();
        }
        string line="";
        while (getline(file, line))
        {
            vector<string> row;
            istringstream is(line);
            for(;;)
            {
                string s = getstr(is);
                row.push_back(s);
                if(is.eof())
                    break;
            }
            res.push_back(row);
        }

        return res;
    }

public:
    unsigned r() const { return (*this).size();}
    unsigned c(unsigned r) const
    {
        assert(r < this->size());
        return (*this)[r].size();
    }
    const string& operator()(unsigned i, unsigned j) const
    {
        ostringstream e;
        if(i>(*this).size())
        {
            e << "csv does not have" << i+1 << "so many lines";
            throw e.str();
        }
        if(j>(*this)[i].size())
        {
            e << "Line " << i+1 << " does not have" << j+1 << "so many columns";
            throw e.str();
        }
        return (*this)[i][j];
    }
    csv(const string& fn) : vector<vector<string>>(getvecs(fn)) {}
};

template<char sep>
inline void csvout(ostream& os, const vector<string>& s)
{
    for(unsigned i=0; i<s.size(); i++)
        csvout<sep>(os, s[i]);
    os << endl;
}

template<typename T, char sep>
inline void csvout(ostream& os, const vector<T>& s)
{
    for(unsigned i=0; i<s.size(); i++)
        os << s[i] << sep;
    os << endl;
}

class matrix
{
public:
    matrix(unsigned r=0, unsigned c=0, double iv=0.0) :
      fx(r,vector<double>(c,iv)) {}

    matrix(const vector<vector<double>>& x): fx(x)
    {
        assert(x.size()>0);
        unsigned m=x[0].size();
        assert(m>0);
        for(unsigned i=1; i<x.size(); i++)
            assert(x[i].size() == m);
    }
    unsigned c() const { return fx[0].size(); }
    unsigned r() const { return fx.size(); }
    double& operator () ( unsigned i, unsigned j )
    { assert(i<r()); assert(j<c()); return fx[i][j]; }
    double operator () ( unsigned i, unsigned j ) const
    { assert(i<r()); assert(j<c()); return fx[i][j]; }

    template <char sep>
    void csvout(ostream& os, bool trans=false, unsigned emptycols=0) const
    {
        if(!trans)
            for(unsigned int i=0; i<r(); i++)
            {
                for(unsigned  j=0; j<emptycols; j++)
                    os << sep;
                for(unsigned  j=0; j<c(); j++)
                    os << fx[i][j] << sep;
                os << endl;
            }
        else
            for(unsigned int j=0; j<c(); j++)
            {
                for(unsigned  j=0; j<emptycols; j++)
                    os << sep;
                for(unsigned int i=0; i<r(); i++)
                    os << fx[i][j] << sep;
                os << endl;
            }
    }
    const vector<double>& operator[] (unsigned i) const
    {
        assert(i<r());
        return fx[i];
    }

protected:
    vector<vector<double>> fx;
};

inline matrix operator * (double a, const matrix& m)
{
    matrix r(m.r(),m.c());
    for(unsigned i=0; i<m.r();i++)
        for(unsigned j=0; j<m.c();j++)
            r(i,j) = a * m(i,j);
    return r;
}

inline matrix operator + (const matrix& a, const matrix& b)
{
    assert(a.r()==b.r());
    assert(a.c()==b.c());
    matrix r(a.r(),a.c());
    for(unsigned i=0; i<a.r();i++)
        for(unsigned j=0; j<a.c();j++)
            r(i,j) = a(i,j) + b(i,j);
    return r;
}

inline double cdist(const matrix& a, const matrix& b)
{
    assert(a.r()==b.r());
    assert(a.c()==b.c());
    double sum = 0;
    for(unsigned i=0;i<a.r();i++)
        for(unsigned j=0;j<a.c();j++)
            sum += fabs(a(i,j)-b(i,j));
    return sum;
}

inline double minoffdiag(const matrix& s)
{
    assert(s.r()==s.c());
    double m=HUGE_VAL;
    for(unsigned i=0; i<s.r(); i++)
        for(unsigned j=0; j<s.r(); j++)
            if(i!=j && s(i,j) < m)
                m = s(i,j);
    return m;
}

inline double maxoffdiag(const matrix& s)
{
    assert(s.r()==s.c());
    double m=-HUGE_VAL;
    for(unsigned i=0; i<s.r(); i++)
        for(unsigned j=0; j<s.r(); j++)
            if(i!=j && s(i,j) > m)
                m = s(i,j);
    return m;
}

inline matrix normoffdiag(const matrix& aw)
{
    double sw = 0.0;
    for(unsigned int i=0; i<aw.r(); i++)
        for(unsigned int j=0; j<aw.c(); j++)
            if(i!=j)
                sw += aw(i,j);

    double avew = sw / (aw.r()*(aw.c()-1));

    return 1.0 / avew * aw;
}

inline double aveoffdiag(const matrix& s,
           const vector<bool>& muted = vector<bool>())
{
    assert(s.r()==s.c());
    assert(muted.size()==0 || muted.size()==s.r());
    unsigned n = muted.size() ?  countfalses(muted) : s.r();
    assert(n > 0);
    if(n=1)
        return 0;
    double sum = 0.0;
    for(unsigned i=0; i<n; i++)
        for(unsigned j=0; j<n; j++)
            if(i!=j && (muted.size()==0 || (!muted[i] && !muted[j])) )
                sum += s(i,j);
    return sum / ( n * (n-1) );
}

inline matrix relchangeoffdiag(const matrix& nm, const matrix& om)
{
    assert(nm.r()==nm.c());
    assert(om.r()==om.c());
    assert(om.r()==nm.c());
    unsigned n=nm.r();
    matrix res(n,n,0);
    for(unsigned i=0; i<n; i++)
        for(unsigned j=0; j<n; j++)
            if(i != j)
                res(i,j) = nm(i,j) / om(i,j);
    return res;
}


} // namespace

#endif // CCTOOLS_HPP
