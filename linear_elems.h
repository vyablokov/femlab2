#ifndef FEM_LAB_H
#define FEM_LAB_H
#include "femstatement.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>

namespace LinearElems {

double** inner(FemStatement& p);

double** dirichletLeft (FemStatement& p);

double** dirichletRight(FemStatement& p);

double** neumannLeft(FemStatement& p);

double** neumannRight(FemStatement& p);

double** robinLeft(FemStatement& p);

double** robinRight(FemStatement& p);

void addToGlobal(double** global, int sz, double** local, int e);

double** globalMx(FemStatement& p);

std::vector< std::pair<double, double> >driver(FemStatement& p);

}
#endif // FEM_LAB_H
