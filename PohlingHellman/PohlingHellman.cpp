//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

/*
 * TODO the g has to have change (check note) , before going to next iteration
 * TODO for each factor => inner loop for i = 0 => exponent ;
 */
NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ g, h, x, divisor;
    std::vector<NTL::ZZ> allXi;
    long y = 0;
    for (const auto &factor:factors) {
        g = NTL::PowerMod(alpha, (N / factor.result), N);
        h = NTL::PowerMod(beta, (N / factor.result), N);
        /*
        * Search for x
        */
        g = NTL::PowerMod(g, NTL::power(factor.prime, factor.exponent - 1), N);
        h = NTL::PowerMod(h, NTL::power(factor.prime, factor.exponent - 1), N);
        PollarRho pollarRho(g, h, N, factor.result);
        x = pollarRho.searchXParallelPollard(); //TODO it takes to long to finish
        allXi.push_back(x);
        for (long i = factor.exponent - 2; i > 0; i++) {
            divisor = 0;
            y = 0;
            std::cout << allXi.size() << " " << allXi.at(0) << std::endl;
            for (auto const &exponent : allXi) {
                divisor = divisor + NTL::MulMod(exponent, NTL::power(factor.prime, y), N);
                y++;
            }
            g = g / divisor;
            h = h / divisor;
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
