//
// Created by mauri on 10.05.18.
//

#include "PohlingHellman.h"
#include "../PollardRho/PollarRho.h"

NTL::ZZ PohlingHellman::searchResult() {
    NTL::ZZ g, h, x, tmpExp, pr;
    NTL::ZZ result = NTL::ZZ(0);
    NTL::ZZ q = N - 1;

    std::vector<NTL::ZZ> allXi;
    std::vector<Congruence> congruences;

    NTL::ZZ tmpG, tmpH;
    for (auto &factor : factors) {
        pr = factor.prime;
//        tmpExp = q / NTL::PowerMod(pr, factors.at(i).exponent, N);
        h = NTL::PowerMod(beta, q / factor.result, N);
        g = NTL::PowerMod(alpha, q / factor.result, N);
        std::cout << h << " <- h g-> " << g << std::endl;
        NTL::ZZ tmp = NTL::ZZ(1);
        NTL::ZZ tmpX = NTL::ZZ(0);
        x = 0;
        tmpG = NTL::PowerMod(g, NTL::power(factor.prime, factor.exponent - 1), N);
        for (long e = 1; e <= factor.exponent; e++) {
            tmpH = NTL::PowerMod(NTL::MulMod(NTL::PowerMod(NTL::InvMod(g, N), x, N), h, N),
                                 NTL::power(factor.prime, factor.exponent - e), N);
            std::cout << tmpG << "<- tmpG tmpH-> " << tmpH << std::endl;
//            PollarRho pollarRho(tmpG, tmpH, N, pr);
//            Pollard pollard(N, pr, tmpG, tmpH);
            tmpX = PollardRhoDLP(N, pr, tmpG, tmpH);
//            tmpX = pollard.findX();
            std::cout << tmpX << " tmpx \n";
//            x += tmpX * NTL::PowerMod(pr, e, N);
            x += tmpX * NTL::power(pr, e - 1);
            allXi.push_back(x);
        }
        congruences.push_back({x, factor.result});
        tmpExp = 0;
        x_factors.push_back(x);
        std::cout << g << "^" << x << "=" << h << " mod" << N << std::endl;
    }
    std::cout << "I have finished to calculate x for each factor\n";

    return SolveCongruences(congruences);
}

NTL::ZZ PohlingHellman::SolveCongruences(std::vector<Congruence> congruences) {
    NTL::ZZ sum = NTL::ZZ(0);
    NTL::ZZ product = NTL::ZZ(1);
    NTL::ZZ gcd, z, k;
    for (auto &congruence : congruences) {
        product *= congruence.result;
    }

    for (auto &congruence : congruences) {
        NTL::ZZ p = product / congruence.result;
        NTL::XGCD(gcd, z, k, p, congruence.result);
//        sum += congruence.x * NTL::InvMod(p % congruence.result, congruence.result) * N;
        sum += z * p * congruence.x;
    }

    return sum % product;
}


NTL::ZZ PohlingHellman::bruteForce(NTL::ZZ base, NTL::ZZ goal, NTL::ZZ up) {
    std::cout << "Parameter " << base << "\n " << goal << "\n " << up << std::endl;
    NTL::ZZ tmpExp;
    for (tmpExp = NTL::ZZ(0); tmpExp <= up; tmpExp++) {
        if (NTL::PowerMod(base, tmpExp, N) == goal) {
            return tmpExp;
        } else if (NTL::PowerMod(base, tmpExp, N) > goal) {
            std::cout << tmpExp << std::endl;
            throw std::invalid_argument("ERR/something horrible happened");
        }
    }
    std::cout << tmpExp << " " << goal - base << std::endl;
    throw std::invalid_argument("ERR/Brute force failed");
}
