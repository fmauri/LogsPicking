//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ g, h, x;
    for (const auto &factor:factors) {
        g = alpha ^ (N / factor.result);
        h = beta ^ (N / factor.result);
        PollarRho pollarRho(g, h, factor.result);
        x = pollarRho.searchXParallelPollard();
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

