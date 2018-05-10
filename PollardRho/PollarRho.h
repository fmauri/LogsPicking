//
// Created by mauri on 10.05.18.
//
#include <NTL/ZZ.h>
#include "../DiscreteLog.h"

#ifndef LOGSPICKING_POLLARRHO_H
#define LOGSPICKING_POLLARRHO_H

struct value {
    NTL::ZZ a, b;
};

class PollarRho : public DiscreteLog {
public:
    PollarRho() {
        length = 45;
        alpha = NTL::ZZ(2);
        x = NTL::RandomLen_ZZ(15);
        beta = NTL::PowerMod(alpha, x, N);
    }

    PollarRho(NTL::ZZ alpha, NTL::ZZ beta, NTL::ZZ N) {
        this->alpha = alpha;
        this->beta = beta;
        this->N = N;
    }

    NTL::ZZ searchXParallelPollard();

private:
    inline void new_xab(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b,
                        const NTL::ZZ &N, const NTL::ZZ &n, const NTL::ZZ &alpha, const NTL::ZZ &beta);

};


#endif //LOGSPICKING_POLLARRHO_H
