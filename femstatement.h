#ifndef FEMSTATEMENT_H
#define FEMSTATEMENT_H
#include <vector>
#include <iostream>

struct Condition {
    enum Kind {dirichlet, neumann, robin};
    Kind kind;
    std::vector<double> vals;
    Condition() {}
};

struct FemStatement
{
    FemStatement(double a, double b, double _m, double _k, double _c, double _f)
        : origin(a), end(b), elemsNum(1), m(_m), k(_k), c(_c), f(_f)
    {
//        std::cout << origin <<" : "<< end << std::endl;
//        std::cout << elemsNum <<' '<< std::endl
//                  << m <<' '<< k <<' '<< c <<' '<< f << std::endl;
    }
    virtual int nodesNum() const {
        return elemsNum+1;
    }
    virtual double elemLen() const {
        return (end - origin)/elemsNum;
    }
    virtual double subElemsLen() const {
        return elemLen();
    }
    double len () const {
        return end - origin;
    }
    double origin;
    double end;
    int elemsNum;
    double m;
    double k;
    double c;
    double f;
    Condition leftCond;
    Condition rightCond;
};

struct FemCubicStatement : public FemStatement {


    // FemStatement interface
public:
    FemCubicStatement(double a, double b, double _m, double _k, double _c, double _f)
        : FemStatement(a, b, _m, _k, _c, _f)
    {}
    int nodesNum() const {
        return 3*elemsNum + 1;
    }
    double subElemsLen() const {
        return elemLen()/3.0;
    }
};


#endif // FEMSTATEMENT_H
