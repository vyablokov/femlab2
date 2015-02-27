#ifndef FEM_COMMON_H
#define FEM_COMMON_H
#include "fstream"

void printMx(double ** mx, int sz, std::ostream &out);

double** allocateMx(int rows, int cols);


#endif // FEM_COMMON_H
