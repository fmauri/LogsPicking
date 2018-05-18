#include "pollard.hpp"
#include <iostream>
#include <NTL/ZZ.h>

using std::cout;
using std::cerr;
using std::endl;
using std::pair;

Pollard::Pollard(ZZ p, ZZ q, ZZ g, ZZ y)
        : p{p}, q{q}, g{g}, y{y} {

}

ZZ Pollard::findX() {
    if (g == ZZ(2)) {
        if (y == g) {
            return ZZ(1);
        } else {
            return ZZ(0);
        }
    }
    if (g == y) {
        return ZZ(1);
    }

    ZZ A = conv<ZZ>(1), B = conv<ZZ>(1);
    size_t i = 0;
    pair<ZZ, ZZ> alpha(NTL::ZZ(0), NTL::ZZ(0));
    pair<ZZ, ZZ> beta(NTL::ZZ(0), NTL::ZZ(0));

    do {
        i++;
        nextStep(A, alpha);
        nextStep(B, beta);
        nextStep(B, beta);
    } while (A != B);

    ZZ diff11 = SubMod(alpha.first, beta.first, q);
    ZZ diff22 = SubMod(beta.second, alpha.second, q);
    return MulMod(diff11, InvMod(diff22, q), q);
}

void Pollard::nextStep(ZZ &a, pair<ZZ, ZZ> &alpha) {
    int modulo = conv<ZZ>(a) % 3;

    if (modulo == 1) {
        a = MulMod(a, y, p);
        alpha.second = AddMod(alpha.second, conv<ZZ>(1), q);
        return;
    }

    if (modulo == 2) {
        a = MulMod(a, a, p);
        alpha.first = MulMod(alpha.first, conv<ZZ>(2), q);
        alpha.second = MulMod(alpha.second, conv<ZZ>(2), q);
        return;
    }

    if (modulo == 0) {
        a = MulMod(a, g, p);
        alpha.first = AddMod(alpha.first, conv<ZZ>(1), q);
        return;
    }
}

