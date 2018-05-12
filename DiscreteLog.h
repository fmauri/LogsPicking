//
// Created by mauri on 10.05.18.
//
#include <NTL/ZZ.h>

#ifndef LOGSPICKING_DISCRETELOG_H
#define LOGSPICKING_DISCRETELOG_H

class DiscreteLog {
protected:
    long length = 40;
    NTL::ZZ n = NTL::GenGermainPrime_ZZ(length, 90); //order of g
    NTL::ZZ N = n * 2 + 1;
    NTL::ZZ alpha;
    NTL::ZZ x;
    NTL::ZZ beta;
public:
    inline void printLog();

};

void DiscreteLog::printLog() {
    std::cout << alpha << "^" << x << " = " << beta << "mod(" << N << ")\n";
}

#endif //LOGSPICKING_DISCRETELOG_H
