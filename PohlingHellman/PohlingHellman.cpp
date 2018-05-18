//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ g, h, x, tmpExp;
    std::vector<NTL::ZZ> allXi;
    h = beta;
    PollarRho pollarRho(alpha, beta, N, n);
    NTL::ZZ q = NTL::ZZ(N - 1);
    for (const auto &factor:factors) {
        for (long e = 1; e < factor.exponent; e++) {
            tmpExp = (N - 1) / NTL::PowerMod(factor.prime, factor.exponent, N);
            h = NTL::PowerMod(h, tmpExp, N);
            g = NTL::PowerMod(alpha, tmpExp, N);
//            x = bruteForce(g, h); //TODO maybe needed to be calculated differently
            pollarRho.setNewValues(g, h, q);
            x = pollarRho.searchXParallelPollard();
            allXi.push_back(x);
            h = NTL::MulMod(beta, NTL::InvMod(NTL::PowerMod(factor.prime, x, N), N), N);
        }
        tmpExp = 0;
        for (auto const &item:allXi) {
            x += item * NTL::PowerMod(alpha, tmpExp, N);
            tmpExp++;
        }
        x_factors.push_back(x);
        std::cout << g << "^" << x << "=" << h << " mod" << N;
    }
    std::cout << "I have finished to calculate x for each factor\n";
//    for (unsigned long i = 0; i < x_factors.size() - 1; i++) {
//        result = result +
//                 NTL::CRT(x_factors.at(i), x_factors.at(i + 1), factors.at(i).result, factors.at(i + 1).result);
//    }
    NTL::ZZ M, y;
    for (unsigned long i = 0; i < x_factors.size() - 1; i++) {
        M = (N - 1) / x_factors.at(i);
        y = NTL::InvMod(M, NTL::ZZ(4));
        result = result + M * y;
    }
    result = result % (N - 1);
    result = result % ((N - 1) / Q);
    return result;
}

NTL::ZZ PohlingHellman::bruteForce(NTL::ZZ num, NTL::ZZ goal) {
    NTL::ZZ tmpExp, tmpResult;
    tmpExp = 1;
    do {
        tmpResult = NTL::PowerMod(num, tmpExp, N);
        ++tmpExp;
    } while (tmpResult != goal);
    return tmpExp;
}
