//
// Created by mauri on 10.05.18.
//

#ifndef LOGSPICKING_POHLINGHELLMAN_H
#define LOGSPICKING_POHLINGHELLMAN_H
#define Factors std::vector<Factor>

#include <vector>
#include "../DiscreteLog.h"

static long factorsLength = 40;

struct Factor {
    NTL::ZZ prime;
    long exponent;
    NTL::ZZ result; //TODO more worth calculate every time or o salve it???????
};

class PohlingHellman : public DiscreteLog {
public:
    PohlingHellman() {
        NTL::ZZ tmpP;
        long exponent;
        NTL::ZZ tmpR;
        do {
            Q = NTL::NextPrime(256, 90);
            factors.clear();
            N = Q;
            exponent = 8 + (rand() % (12 - 8 + 1));
            N = N * 2 ^ exponent;
            for (int i = 0; i < 8; i++) {
                tmpP = NTL::NextPrime(factorsLength, 90);
                exponent = 8 + (rand() % (12 - 8 + 1));
                tmpR = tmpP ^ exponent;
                Factor tmpF{tmpP, exponent, tmpR};
                factors.push_back(tmpF);
                N = N * (tmpR);
            }
            N++; // +1 to make it odd and this should be prime
        } while (NTL::ProbPrime(N, 1000));
        alpha = NTL::RandomBnd(N); // up to N-1
        x = NTL::RandomBnd(N - 1) + 1;
        beta = NTL::PowerMod(alpha, x, N);
    }

    NTL::ZZ searchResult();

private:
    NTL::ZZ Q;
    Factors factors;
};

#endif //LOGSPICKING_POHLINGHELLMAN_H
