//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ g, h, x, tmpExp;
    std::vector<NTL::ZZ> allXi;
    long y = 0;
    for (const auto &factor:factors) {
        tmpExp = (N - 1) / factor.result;
        g = NTL::PowerMod(alpha, tmpExp, N);
        h = NTL::PowerMod(beta, tmpExp, N);
        /*
        * Search for x
        */
        g = NTL::PowerMod(g, NTL::power(factor.prime, factor.exponent - 1), N);
        h = NTL::PowerMod(h, NTL::power(factor.prime, factor.exponent - 1), N);
        PollarRho pollarRho(g, h, N, N - 1);
        std::cout << std::endl;
        std::cout << g << "^x=" << h << "mod" << N;
        x = pollarRho.searchXParallelPollard();
        allXi.push_back(x);
#pragma omp parallel for schedule(dynamic)
        for (long i = factor.exponent - 2; i > 0; i++) {
            tmpExp = 0;
            y = 0;
            for (auto const &exponent : allXi) {
                tmpExp = tmpExp + NTL::MulMod(exponent, NTL::power(factor.prime, y), N);
                y++;
            }
            h = NTL::PowerMod(NTL::MulMod(h, NTL::PowerMod(h, -tmpExp, N), N), i, N);
            pollarRho.setNewValues(g, h, factor.result);
            x = pollarRho.searchXParallelPollard();
            allXi.push_back(x);
        }
        x = 0;
        for (auto const &item:allXi) {
            x += item;
        }
        x_factors.push_back(x);
        std::cout << "Factor done\n";
    }
    std::cout << "I have finished to calculate x for each factor\n";
    for (unsigned long i = 0; i < x_factors.size() - 1; i++) {
        result = result +
                 NTL::CRT(x_factors.at(i), x_factors.at(i + 1), factors.at(i).result, factors.at(i + 1).result);
    }
//    result = calcCRT();
    result = result % ((N - 1) / Q);
    return result;
}
