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
        g = alpha ^ (N / factor.result);
        h = beta ^ (N / factor.result);
        /*
         * Search for x
         */
        g = NTL::PowerMod(g, NTL::power(factor.prime, factor.exponent - 1), N);
        h = NTL::PowerMod(h, NTL::power(factor.prime, factor.exponent - 1), N);
        PollarRho pollarRho(g, h, N, factor.prime);
        x = pollarRho.searchXParallelPollard();
        allXi.push_back(x);
        for (long i = factor.exponent - 2; i > 0; i++) {
            divisor = 0;
            y = 0;
            for (auto const &exponent : allXi) {
                divisor = divisor + (exponent * NTL::power(factor.prime, y));
                y++;
            }
            g = g / divisor;
            h = h / divisor;
            pollarRho.setNewValues(g, h);
            x = pollarRho.searchXParallelPollard();
            allXi.push_back(x);
        }
        x = 0;
        for (auto const &item:allXi) {
            x += item;
        }
        x_factors.push_back(x);
    }
    result = calcCRT();
    return result;
}

NTL::ZZ PohlingHellman::calcCRT() {
    NTL::ZZ result = NTL::ZZ(1); // Initialize result

    // As per the Chinese remainder theorem,
    // this loop will always break.
    while (true) {
        // Check if remainder of result % num[j] is
        // rem[j] or not (for all j from 0 to k-1)
        unsigned long i;
        for (i = 0; i < x_factors.size(); i++) {
            if (result % x_factors.at(i) != factors.at(i).result) {
                break;
            }
        }
        // If all remainders matched, we found result
        if (i == x_factors.size()) {
            break;
//            return result;
        }

        // Else try next number
        result++;
    }

    return result;
}
