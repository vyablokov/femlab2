#ifndef CUBIC_ELEMS_H
#define CUBIC_ELEMS_H

#include "femstatement.h"
#include <vector>
#include <utility>


namespace CubicElems {

double** inner(FemCubicStatement &p);

double** dirichletLeft (FemCubicStatement& p);

double** dirichletRight(FemCubicStatement& p);

double** neumannLeft(FemCubicStatement& p);

double** neumannRight(FemCubicStatement& p);

double** robinLeft(FemCubicStatement& p);

double** robinRight(FemCubicStatement& p);

void addToGlobal(double** global, int sz, double** local, int e);

double** globalMx(FemCubicStatement& p);

std::vector<std::pair<double, double> > driver(FemCubicStatement& p);

}

#endif // CUBIC_ELEMS_H
