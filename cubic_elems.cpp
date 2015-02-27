#include "cubic_elems.h"
#include <fem_common.h>
#include "femstatement.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "slae.h"
#include <vector>
#include <utility>

namespace CubicElems {

double** inner(FemCubicStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(4,5);

    locMx[0][0] = 8.0*c*l/105.0 - 37.0*m/(10.0*l) - 0.5*k;
    locMx[0][1] = 57*k/80.0 + 189.0*m/(40.0*l) + 33.0*c*l/560.0;
    locMx[0][2] = -0.3*k - 27.0*m/(20.0*l) - 3*c*l/140.0;
    locMx[0][3] = 7.0*k/80.0 + 13.0*m/(40.0*l) + 19.0*c*l/1680.0;
    locMx[0][4] = f*l/8.0;

    locMx[1][0] = 189.0*m/(40.0*l) - 57.0*k/80.0 + 33.0*c*l/560.0;
    locMx[1][1] = 27*(c*l*l - 28.0*m)/(70*l);
    locMx[1][2] = 81.0*k/80.0 + 297.0*m/(40.0*l) - 27*c*l/560.0;
    locMx[1][3] = -0.3*k - 27*m/(20.0*l) - 3*c*l/140.0;
    locMx[1][4] = 3*f*l/8.0;

    locMx[2][0] = 0.3*k - 27.0*m/(20.0*l) - 3*c*l/140.0;
    locMx[2][1] = 297.0*m/(40.0*l) - 81.0*k/80.0 - 27.0*c*l/560.0;
    locMx[2][2] = 27.0*(c*l*l - 28.0*m)/(70.0*l);
    locMx[2][3] = 57*k/80.0 + 189*m/(40.0*l) + 33*c*l/560.0;
    locMx[2][4] = 3*f*l/8.0;


    locMx[3][0] = 13*m/(40.0*l) - 7.0*k/80.0 + 19.0*c*l/1680.0;
    locMx[3][1] = 0.3*k - 27.0*m/(20.0*l) - 3*c*l/140.0;
    locMx[3][2] = 189.0*m/(40.0*l) - 57.0*k/80.0 + 33.0*c*l/560.0;
    locMx[3][3] = 0.5*k - 37.0*m/(10.0*l) + 8*c*l/105.0;
    locMx[3][4] = f*l/8.0;

    return locMx;
}

double** dirichletLeft (FemCubicStatement &p)
{
    double** locMx{nullptr};
    locMx = CubicElems::inner(p); // get original inner element local matrix

    // toss first ([0]) column to the right side multyplied by u0 = u(x0)
    for (int i = 0; i < 4; i++) {
        locMx[i][4] -= locMx[i][0]*p.leftCond.vals[0];
        locMx[i][0] = 0.0;
    }
    locMx[0][0] = -p.m;

    return locMx;
}

double** dirichletRight(FemCubicStatement &p)
{
    double** locMx{nullptr};
    locMx = CubicElems::inner(p);

    // toss fourth ([3]) column to the right side multyplied by u1 = u(x1)
    for (int i = 0; i < 3; i++) {
        locMx[i][4] -= locMx[i][3]*p.rightCond.vals[0];
        locMx[i][3] = 0.0;
    }
    locMx[3][3] = p.m;

    return locMx;
}

double** neumannLeft(FemCubicStatement &p)
{
    double** locMx{nullptr};
    locMx = CubicElems::inner(p);

    // add {m * u'(x0); 0...0} to the right side
    locMx[0][4] += p.m * p.leftCond.vals[0];

    return locMx;
}

double** neumannRight(FemCubicStatement &p)
{
    double** locMx{nullptr};
    locMx = CubicElems::inner(p);

    // add {0...0; -m * u'(x1)} to the right side
    locMx[3][4] -= p.m * p.rightCond.vals[0];

    return locMx;
}

double** robinLeft(FemCubicStatement &p)
{
    double k1{p.leftCond.vals[0] / p.leftCond.vals[1]}; // alpha / beta
    double k2{p.leftCond.vals[2] / p.leftCond.vals[1]}; // gamma / beta
    double** locMx{nullptr};

    locMx = CubicElems::inner(p);

    locMx[0][0] += p.m * k1;
    locMx[0][4] += p.m * k2;

    return locMx;
}

double** robinRight(FemCubicStatement &p)
{
    double k1{p.leftCond.vals[0] / p.leftCond.vals[1]};
    double k2{p.leftCond.vals[2] / p.leftCond.vals[1]};
    double** locMx{nullptr};

    locMx = CubicElems::inner(p);
    locMx[3][3] -= p.m * k1;
    locMx[3][4] -= p.m * k2;

    return locMx;
}

void addToGlobal(double** global, int sz, double** local, int e)
{
    int d{3};

    for (int i = 0; i < d+1; i++) {
        for (int j = 0; j < d+1; j++)
            global[d*e+i][d*e+j] += local[i][j];
        global[d*e+i][sz] += local[i][d+1];
    }
}

double** globalMx(FemCubicStatement &p)
{
    double **left, **right, **in;

    switch (p.leftCond.kind)
    {
    case Condition::Kind::dirichlet:
        left = dirichletLeft(p);
        break;
    case Condition::Kind::neumann:
        left = neumannLeft(p);
        break;
    case Condition::Kind::robin:
        left = robinLeft(p);
    default:
        break;
    }

    switch (p.rightCond.kind)
    {
    case Condition::Kind::dirichlet:
        right = dirichletRight(p);
        break;
    case Condition::Kind::neumann:
        right = neumannRight(p);
        break;
    case Condition::Kind::robin:
        right = robinRight(p);
    default:
        break;
    }

    in = inner(p);

    double** mx;
    mx = allocateMx(p.nodesNum(), p.nodesNum()+1);

    addToGlobal(mx, p.nodesNum(), left, 0);

    for (int e = 1; e < p.elemsNum-1; e++)
        addToGlobal(mx, p.nodesNum(), in, e);
    addToGlobal(mx, p.nodesNum(), right, p.elemsNum-1);

    std::ofstream file("file");
    printMx(mx, p.nodesNum(), file);
    return mx;
}

std::vector<std::pair<double, double> > driver(FemCubicStatement &p)
{
    using namespace std;
    SLAE s;
    s.insertMatrix(globalMx(p), p.nodesNum());
    s.solve();
    double *S = new double[p.nodesNum()];

    s.copySoln(S);
    // simple postprocessing
    if (p.leftCond.kind == Condition::Kind::dirichlet)
        S[0] = p.leftCond.vals[0];

    if (p.rightCond.kind == Condition::Kind::dirichlet)
        S[p.nodesNum()-1] = p.rightCond.vals[0];

    // building a QVector<QPointF> to make samples for QwtPlotCurve
    vector< pair<double, double> > vec;
    for (int i = 0; i < p.nodesNum(); i++)
        vec.push_back( pair<double, double>(p.origin + p.subElemsLen()*i, S[i]) );

    delete [] S;
    return vec;
}

}

