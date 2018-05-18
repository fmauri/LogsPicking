//
// Created by mauri on 10.05.18.
//

#ifndef LOGSPICKING_POLLARRHO_H
#define LOGSPICKING_POLLARRHO_H

#include <NTL/ZZ.h>
#include "../DiscreteLog.h"

struct value {
    NTL::ZZ a, b;
};

class PollardRho : public DiscreteLog {
public:
    PollardRho() {
        length = 45;
        alpha = NTL::ZZ(2);
        x = NTL::RandomLen_ZZ(15);
        beta = NTL::PowerMod(alpha, x, N);
    }

    PollardRho(NTL::ZZ alpha, NTL::ZZ beta, NTL::ZZ N, NTL::ZZ n) {
        this->alpha = alpha;
        this->beta = beta;
        this->N = N;
        this->n = n;
    }

    NTL::ZZ setNewValues(NTL::ZZ alpha, NTL::ZZ beta, NTL::ZZ ord) {
        this->alpha = alpha;
        this->beta = beta;
        this->n = ord;
    };

    NTL::ZZ searchXParallelPollard();

private:
    inline void new_xab(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b,
                        const NTL::ZZ &N, const NTL::ZZ &n, const NTL::ZZ &alpha, const NTL::ZZ &beta);

};


#endif //LOGSPICKING_POLLARRHO_H
