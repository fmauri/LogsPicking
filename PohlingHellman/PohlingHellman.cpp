//
// Created by mauri on 02.06.18.
//

#include "PohlingHellman.h"
#include "pollard_rho_dlp.hpp"

void PohlingHellman::GenerateChallenge() {
    while (true) {
        long k = NTL::RandomBnd(5) + 8;
        double powerBrowsingLimit = pow(5, k);

        Q = NTL::GenPrime_ZZ(NTL::RandomBnd(32) + 256) * NTL::GenPrime_ZZ(NTL::RandomBnd(32) + 256);
        primes.push_back(NTL::ZZ(2));

        for (int i = 1; i < k; ++i) {
            primes.push_back(NTL::GenGermainPrime_ZZ(NTL::RandomBnd(11) + 40));
        }

        for (int i = 0; i < powerBrowsingLimit; ++i) {
            int browsePoint = i;

            exponents.clear();
            P = Q;

            for (int j = 0; j < k; ++j) {
                exponents.push_back(browsePoint % 5 + 3);
                P *= NTL::power(primes[j], exponents[j]);
                browsePoint /= 5;
            }

            N = P;

            if (++P % 2 == 1 && !NTL::MillerWitness(P, NTL::ZZ(32))) {
                alpha = FindRandomElementOfGroup(P);
                x = FindRandomElementOfGroup(P);
                beta = NTL::PowerMod(alpha, x, P);

                return;
            }
        }
    }
}

NTL::ZZ PohlingHellman::calculateReduction() {
    std::vector<CongruentPoint> congruences;
    NTL::ZZ x, q, e, a, b, c, l;

    // Find the congruence using exhaustive method for p0^e0. (p0 = 2)
    NTL::ZZ qe2 = NTL::PowerMod(primes[0], exponents[0], P);
    NTL::ZZ g2 = NTL::PowerMod(alpha, N / qe2, P);
    NTL::ZZ y2 = NTL::PowerMod(beta, N / qe2, P);
    NTL::ZZ x2 = ExhaustiveDLP(P, qe2, g2, y2);
    congruences.push_back({x2, qe2});

    // Find the congruences using Pollard-Rho method for p1^e1, ..., pk^ek.
    for (int i = 1; i < primes.size(); ++i) {
        x = 0;
        q = primes[i];
        e = exponents[i];

        c = 1;
        l = 0;

        a = NTL::PowerMod(alpha, N / q, P);

        for (int j = 0; j < e; ++j) {
            c = j > 0 ? NTL::MulMod(c, NTL::PowerMod(alpha, NTL::MulMod(l, NTL::PowerMod(q, j - 1, P), P), P), P)
                      : NTL::ZZ(1);
            b = NTL::PowerMod(NTL::MulMod(beta, NTL::InvMod(c, P), P), N / NTL::PowerMod(q, j + 1, P),
                              P);
            l = PollardRhoDLP(P, q, a, b);

            x += l * NTL::PowerMod(q, j, P);
        }

        congruences.push_back({x, NTL::PowerMod(q, e, P)});
    }

    return SolveCongruences(congruences);
}

NTL::ZZ PohlingHellman::ExhaustiveDLP(NTL::ZZ p, NTL::ZZ q, NTL::ZZ g, NTL::ZZ y) {
    for (int i = 0; i < q; ++i) {
        if (NTL::PowerMod(g, i, p) == y) {
            return NTL::ZZ(i);
        }
    }

    throw std::invalid_argument("Something went horribly wrong.");
}

NTL::ZZ PohlingHellman::SolveCongruences(std::vector<CongruentPoint> congruences) {
    NTL::ZZ sum = NTL::ZZ(0);
    NTL::ZZ product = NTL::ZZ(1);

    for (int i = 0; i < congruences.size(); ++i) {
        product *= congruences[i].P;
    }

    for (int i = 0; i < congruences.size(); ++i) {
        NTL::ZZ p = product / congruences[i].P;
        sum += congruences[i].X * NTL::InvMod(p % congruences[i].P, congruences[i].P) * p;
    }

    return sum % product;
}
