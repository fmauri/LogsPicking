//
// Created by mauri on 02.06.18.
//

#ifndef LOGSPICKING_POHLINGHELLMAN_H
#define LOGSPICKING_POHLINGHELLMAN_H


#include <NTL/ZZ.h>
#include <vector>

struct CongruentPoint {
    NTL::ZZ X;
    NTL::ZZ P;
};

class PohlingHellman {
public:
    void GenerateChallenge();

    NTL::ZZ calculateReduction();

    NTL::ZZ getGoal() { return x % (N / Q); };

    NTL::ZZ getAlphaQ() { return NTL::PowerMod(alpha, Q, P); }

    NTL::ZZ getBetaQ() { return NTL::PowerMod(beta, Q, P); }

    NTL::ZZ getN() { return N; }

private:
    std::vector<NTL::ZZ> primes;
    std::vector<int> exponents;
    NTL::ZZ alpha, x, beta;
    NTL::ZZ Q, N, P;

    inline NTL::ZZ FindRandomElementOfGroup(NTL::ZZ p) { return NTL::RandomBnd(p - 1) + NTL::ZZ(1); }

    NTL::ZZ ExhaustiveDLP(NTL::ZZ p, NTL::ZZ q, NTL::ZZ g, NTL::ZZ y);

    NTL::ZZ SolveCongruences(std::vector<CongruentPoint> congruences);

};


#endif //LOGSPICKING_POHLINGHELLMAN_H
