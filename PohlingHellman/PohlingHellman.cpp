//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ g, h, x;
    std::vector<NTL::ZZ> x_factors;
    for (const auto &factor:factors) {
        g = alpha ^ (N / factor.result);
        h = beta ^ (N / factor.result);
        PollarRho pollarRho(g, h, N);
        x = pollarRho.searchXParallelPollard();
        x_factors.push_back(x);
    }

    return result;
}
