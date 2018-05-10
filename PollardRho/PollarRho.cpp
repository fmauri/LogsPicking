//
// Created by mauri on 10.05.18.
//

#include <omp.h>
#include <map>
#include "PollarRho.h"

void PollarRho::new_xab(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b, const NTL::ZZ &N, const NTL::ZZ &n, const NTL::ZZ &alpha,
                        const NTL::ZZ &beta) {
    switch (x % 3) {
        case 0:
            x = PowerMod(x, 2, N);
            a = MulMod(a, 2, n);
            b = MulMod(b, 2, n);
            break;
        case 1:
            x = MulMod(x, alpha, N);
            a = AddMod(a, 1, n);
            break;
        case 2:
            x = MulMod(x, beta, N);
            b = AddMod(b, 1, n);
            break;
        default:
            break;
    }
}

NTL::ZZ PollarRho::searchXParallelPollard() {
    const NTL::ZZ n = (N - 1) / 2;
    std::map<NTL::ZZ, value> collisions;

    NTL::ZZ x, a, b;
    NTL::ZZ tmp1, tmp2;
    NTL::ZZ result = NTL::ZZ(0);
    int j;
    NTL::ZZ steps = NTL::SqrRoot(M_PI * N / 2) / 8;
    omp_set_num_threads(8);
#pragma omp parallel for private(x, a, b, tmp1, tmp2) shared(result, collisions) schedule(dynamic)
    for (j = 0; j < 1000; j++) {
        if (result == 0) {
            a = NTL::RandomBnd(N);
            b = NTL::RandomBnd(N);
            x = NTL::MulMod(NTL::PowerMod(alpha, a, N), NTL::PowerMod(beta, b, N), N);
            for (NTL::ZZ i = NTL::ZZ(0); i < steps; i++) {
                new_xab(x, a, b, N, n, alpha, beta);
                if ((x % 0xfff) == 0) {
                    if (collisions.find(x) != collisions.end()) {
                        if (result == 0) {
                            NTL::ZZ nominator = NTL::SubMod(a, collisions.at(x).a, n);
                            NTL::ZZ denominator = NTL::SubMod(collisions.at(x).b, b, n);
                            result = NTL::MulMod(NTL::InvMod(denominator, n), nominator, n);
                        }
                    } else {
                        value tmp_result = {a, b};
                        collisions.insert({x, tmp_result});
                    }
                    break;
                }
            }
        }
    }
    return result;
}
