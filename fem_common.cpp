#include "fem_common.h"

#include <iostream>
#include <fstream>
#include <iomanip>

void printMx(double ** mx, int sz, std::ostream& out)
{
    out << std::fixed;
    out << std::setprecision(2);
    for (int i = 0; i < sz; i++) {
        for (int j = 0; j < sz+1; j++)
            out << std::setw(8) << mx[i][j] <<' ';
        out << std::endl;
    }
}

double** allocateMx(int rows, int cols)
{
    double **m{nullptr};

    m = new double* [rows];

    for (int i = 0; i < rows; i++) {
        m[i] = new double [cols];
        for (int j = 0; j < cols; j++)
            m[i][j] = 0.0;
    }

    return m;
}
