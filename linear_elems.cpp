#include "linear_elems.h"
#include <fem_common.h>
#include "femstatement.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "slae.h"
#include <qvector.h>
#include <QPointF>

namespace LinearElems {

double** inner(FemStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5;

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0;
    locMx[1][2] = f*l*0.5;

    return locMx;
}

double** dirichletLeft (FemStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(2, 3);

    locMx[0][0] = -m;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5 - p.leftCond.vals[0]*(-m/l - 0.5*k + c*l/3.0);

    locMx[1][0] = 0.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0;
    locMx[1][2] = f*l*0.5 - p.leftCond.vals[0]*(m/l - 0.5*k + c*l/6.0);


    return locMx;
}

double** dirichletRight(FemStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0;
    locMx[0][1] = 0.0;
    locMx[0][2] = f*l*0.5 - p.rightCond.vals[0]*(m/l + 0.5*k + c*l/6.0);

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = m;
    locMx[1][2] = f*l*0.5 - p.rightCond.vals[0]*(-m/l + 0.5*k + c*l/3.0);

    return locMx;
}

double** neumannLeft(FemStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5 + m*p.leftCond.vals[0];

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0;
    locMx[1][2] = f*l*0.5;

    return locMx;
}

double** neumannRight(FemStatement& p)
{
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};
    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5;

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0;
    locMx[1][2] = f*l*0.5 - m*p.rightCond.vals[0];

    return locMx;
}

double** robinLeft(FemStatement& p)
{
    double k1{p.leftCond.vals[0] / p.leftCond.vals[1]};
    double k2{p.leftCond.vals[2] / p.leftCond.vals[1]};
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};

    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0 + m*k1;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5 + m*k2;

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0;
    locMx[1][2] = f*l*0.5;

    return locMx;
}

double** robinRight(FemStatement& p)
{
    double k1{p.leftCond.vals[0] / p.leftCond.vals[1]};
    double k2{p.leftCond.vals[2] / p.leftCond.vals[1]};
    double m{p.m}, k{p.k}, c{p.c}, f{p.f}, l{p.elemLen()};
    double** locMx{nullptr};

    locMx = allocateMx(2, 3);

    locMx[0][0] = -m/l - 0.5*k + c*l/3.0;
    locMx[0][1] = m/l + 0.5*k + c*l/6.0;
    locMx[0][2] = f*l*0.5;

    locMx[1][0] = m/l - 0.5*k + c*l/6.0;
    locMx[1][1] = -m/l + 0.5*k + c*l/3.0 - m*k1;
    locMx[1][2] = f*l*0.5 - m*k2;

    return locMx;
}

void addToGlobal(double** global, int sz, double** local, int e)
{
    int d{1};

    for (int i = 0; i < d+1; i++) {
        for (int j = 0; j < d+1; j++)
            global[d*e+i][d*e+j] += local[i][j];
        global[d*e+i][sz] += local[i][d+1];
    }
}

double** globalMx(FemStatement& p)
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

    return mx;
}

std::vector<std::pair<double, double> > driver(FemStatement& p)
{
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
    std::vector<std::pair<double, double> > vec;
    for (int i = 0; i < p.nodesNum(); i++)
        vec.push_back(std::pair<double, double>(p.origin + p.elemLen()*i, S[i]) );

    delete [] S;
    return vec;
}

}
