#ifndef SLAE_H
#define SLAE_H
#include <cstdlib>
#include <iostream>
#include <iomanip>

class SLAE
{
    unsigned sz;
    double** mx; // matrix sz*(sz+1)
    double EPSYLON;
    void gauss_fwd (double**, int, int, int);
    void gauss_bwd (double**, int, int, int);
    friend std::ostream& operator << (std::ostream& out, SLAE& s);
public:
    SLAE();
    void solve(double eps = 10.0E-8);
    unsigned size() const {
        return sz;
    }
    void copySoln(double* S)
    {
        if (S == nullptr) exit(1);
        for (unsigned i = 0; i < sz; i++)
            S[i] = mx[i][sz];
    }
    void insertMatrix(double** M, int s)
    {
        sz = s;
        mx = M;
    }

    void insertFree(double* B)
    {
        if (B == nullptr) exit(2);
        for (unsigned i = 0; i < sz; i++)
            mx[i][sz] = B[i];
    }
};



#endif // SLAE_H
