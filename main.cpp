#include "femstatement.h"
#include "linear_elems.h"
#include "cubic_elems.h"
#include <fstream>

int main(/*int argc, char *argv[]*/)
{
    using namespace std;
		// FemCubicStatement constructor params:
		// xbegin,  xend, u'' coef, u' coef, u coef, right side  
    FemCubicStatement p(1.0, 8.0, 11.0, 3.0, 0.0, 7.0);
    p.elemsNum = 20; // elemsNum)))
		 // boundary conditions of the first kind - dirichlet
		 // boundary conditions of the second kind - neumann
    p.leftCond.kind = Condition::Kind::dirichlet; // first kind
    p.leftCond.vals.push_back(-6.0); // u(xbegin) = u0 - push_back(u0)
    p.rightCond.kind  = Condition::Kind::neumann;  // second kind
    p.rightCond.vals.push_back(-6.0); // u'(xend) = u'1 - push_back(u'1)

    auto vec = CubicElems::driver(p);

    ofstream file("results.svd");
    for (auto & elem : vec)
        file << elem.first << ' ' << elem.second << endl;

    return 0;
}
