#include "slae.h"
#include <math.h>

SLAE::SLAE() : EPSYLON(10.0E-8)
{
}

void SLAE::solve(double eps)
{
    EPSYLON = eps;
    gauss_fwd(mx, sz, 0, sz-1);
    gauss_bwd(mx, sz, sz-1, 0);
}

void SLAE::gauss_fwd(double **M, int size, int s, int f)
{
    int nRows = size;
    int nCols = nRows+1;

    /* forward Gauss method step */
    int col, row, keyRow;
    double k;
    for (keyRow = s; keyRow <= f; keyRow++) {
        for (row = keyRow+1; row < nRows; row++) {
            if (fabs(M[row][keyRow]) > EPSYLON) {
                k = -M[row][keyRow] / M[keyRow][keyRow];
                for (col = 0; col < nCols; col++)
                    if ( fabs(M[keyRow][col]) > EPSYLON )
                        M[row][col] += M[keyRow][col] * k;
            }
        }
    }
    return;
}

void SLAE::gauss_bwd(double **M, int size, int s, int f)
{
    int keyRow, nRows = size, nCols = nRows+1, col/*, row*/;
    double k;

    /* backward Gauss method step */
    for (keyRow = s; keyRow >= f; keyRow--) {
        k = 0.0;
        for (col = keyRow+1/*, row = nRows-1*/; col < nCols-1; col++) {
            k += M[keyRow][col] * M[col][nCols-1];
            /*	printf("k += %f * %f ", M[keyRow][col], M[col][nCols-1]);*/
        }
        /*	printf("keyRow = %d k = %f\n", keyRow, k);*/
        M[keyRow][nCols-1] -= k;
        M[keyRow][nCols-1] /= M[keyRow][keyRow];
    }

    return;
}

std::ostream& operator << (std::ostream& out, SLAE& s) {

    out << std::fixed;
    out << std::setprecision(1);
    for (unsigned i = 0; i < s.sz; i++) {
        for (unsigned j = 0; j < s.sz+1; j++)
            out << std::setw(6) << s.mx[i][j] <<' ';
        out << std::endl;
    }
    return out;
}
