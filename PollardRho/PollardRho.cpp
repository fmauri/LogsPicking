//
// Created by mauri on 10.05.18.
//

#include <omp.h>
#include <map>
#include "PollardRho.h"

void PollardRho::new_xab(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b, const NTL::ZZ &N, const NTL::ZZ &n, const NTL::ZZ &alpha,
                         const NTL::ZZ &beta) {
    switch (x % 3) {
        case 0:
//            x = NTL::PowerMod(x, 2, N);
            x = NTL::MulMod(x, x, N);
            a = NTL::MulMod(a, 2, n);
            b = NTL::MulMod(b, 2, n);
            break;
        case 1:
            x = NTL::MulMod(x, alpha, N);
            a = NTL::AddMod(a, 1, n);
            break;
        case 2:
            x = NTL::MulMod(x, beta, N);
            b = NTL::AddMod(b, 1, n);
            break;
        default:
            break;
    }
}

NTL::ZZ PollardRho::searchXParallelPollard() {
    std::map<NTL::ZZ, value> collisions;

    NTL::ZZ x, a, b;
    NTL::ZZ tmp1, tmp2;
    NTL::ZZ result = NTL::ZZ(0);
    int j;
    NTL::ZZ steps = NTL::SqrRoot(M_PI * N / 2) / 8;
    omp_set_num_threads(8);
    if (alpha == 2) {
        if (beta == alpha) {
            return NTL::ZZ(1);
        }/* else {
            return NTL::ZZ(0);
        }*/
    }
    if (alpha == beta) {
        return NTL::ZZ(1);
    }
    omp_set_num_threads(8);
#pragma omp parallel for private(x, a, b, tmp1, tmp2, j) shared(result, collisions) schedule(dynamic)
    for (j = 0; j < omp_get_thread_num(); j++) {
        while (result == 0) {
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
