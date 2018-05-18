//
// Created by mauri on 10.05.18.
//

#ifndef LOGSPICKING_POHLINGHELLMAN_H
#define LOGSPICKING_POHLINGHELLMAN_H
#define Factors std::vector<Factor>

#include <vector>
#include "../DiscreteLog.h"
#include "../PollardRho/PollardRho.h"

static long factorsLength = 35;

struct Factor {
    NTL::ZZ prime;
    long exponent;
    NTL::ZZ result;
};

struct Congruence {
    NTL::ZZ x;
    NTL::ZZ result;
};

class PohlingHellman : public DiscreteLog {
public:
    PohlingHellman() {
        NTL::ZZ tmpP;
        long exponent;
        NTL::ZZ tmpR;
        do {
            long factorsAmount = NTL::RandomBnd(5) + 8;
            factors.clear();
            Q = NTL::GenPrime_ZZ(256);
            N = Q;
            exponent = NTL::RandomBnd(5) + 8;
            tmpP = 2;
            tmpR = NTL::power(tmpP, exponent);
            Factor tmpF1{tmpP, exponent, tmpR};
            factors.push_back(tmpF1);
            N = N * tmpR;
            for (int i = 0; i < factorsAmount; i++) {
                tmpP = NTL::GenGermainPrime_ZZ(factorsLength);
                exponent = NTL::RandomBnd(10) + 3;
                tmpR = NTL::power(tmpP, exponent);
                Factor tmpF{tmpP, exponent, tmpR};
                factors.push_back(tmpF);
                N = N * (tmpR);
            }
            N = N + 1; // +1 to make it odd and this should be prime
        } while (NTL::ProbPrime(N, 1000000));
        alpha = NTL::RandomBnd(N - 2) + 2;
        x = NTL::RandomBnd(N - 2) + 1;
        beta = NTL::PowerMod(alpha, x, N);
//        N = 8101;
//        alpha = 6;
//        x = 6689;
//        beta = 7531;
//        factors.clear();
//        Factor factor{NTL::ZZ(2), 2, NTL::power(NTL::ZZ(2), 2)};
//        Factor factor1{NTL::ZZ(3), 4, NTL::power(NTL::ZZ(3), 4)};
//        Factor factor2{NTL::ZZ(5), 2, NTL::power(NTL::ZZ(5), 2)};
//        factors.push_back(factor);
//        factors.push_back(factor1);
//        factors.push_back(factor2);
    }

    NTL::ZZ searchResult();

    NTL::ZZ bruteForce(NTL::ZZ base, NTL::ZZ goal, NTL::ZZ up);

    NTL::ZZ SolveCongruences(std::vector<Congruence> congruences);

    void printFactors();

private:
    NTL::ZZ Q;
    Factors factors;
    std::vector<NTL::ZZ> x_factors;
};

#endif //LOGSPICKING_POHLINGHELLMAN_H
