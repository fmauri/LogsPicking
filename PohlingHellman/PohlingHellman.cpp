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
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ p, i, prod = NTL::ZZ(1), sum = NTL::ZZ(0);

    for (const auto &factor:factors) {
        prod *= factor.result;
    }

    for (unsigned long y = 0; y < factors.size(); y++) {
        p = prod / factors.at(y).result;
        sum += x_factors.at(y) * mul_inv(p, factors.at(y).result) * p;
    }

    return result % prod;
}

NTL::ZZ PohlingHellman::mul_inv(NTL::ZZ a, NTL::ZZ b) {
    NTL::ZZ b0 = b, t, q;
    NTL::ZZ x0 = NTL::ZZ(0), x1 = NTL::ZZ(1);
    if (b == 1) return NTL::ZZ(1);
    while (a > 1) {
        q = a / b;
        t = b, b = a % b, a = t;
        t = x0, x0 = x1 - q * x0, x1 = t;
    }
    if (x1 < 0) x1 += b0;
    return x1;
}

