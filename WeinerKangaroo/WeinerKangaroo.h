//
// Created by mauri on 10.05.18.
//

#ifndef LOGSPICKING_WEINERKANGAROO_H
#define LOGSPICKING_WEINERKANGAROO_H
#define jumps_info std::vector<NTL::ZZ>

#include <NTL/ZZ.h>
#include <vector>
#include "../DiscreteLog.h"

class WeinerKangaroo : public DiscreteLog {
public:
    WeinerKangaroo() {
        x = NTL::RandomBnd(lower - upper) + upper;
        alpha = setAlpha();
        beta = NTL::PowerMod(alpha, x, N);
        fillSets();
    }

    NTL::ZZ searchCollisions();

private:

    const long length = 45;
    const long numProcess = 8;
    const long numWild = numProcess / 2 + 1; // v
    const long numTamed = numProcess / 2 - 1; // u
    const long GAMMA = 65536;
    const long OMEGA = 0xffff;
    const NTL::ZZ lower = NTL::RandomBnd(N - (n * M_PI / numProcess / 2)); // a
    const NTL::ZZ upper = lower + n * M_PI / numProcess / 2; // b
    const NTL::ZZ median = (upper - lower) / 2; // (a+b)/2
    const unsigned long k = std::floor(NTL::log(upper - lower));
    jumps_info distanceS; // S
    jumps_info jumpSetR; // R
    std::hash<std::string> hash_fn;

    void performStep(NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps);

    std::string zToString(const NTL::ZZ &z);

    inline void fillSets();

    inline void initKangaroo(int ind, NTL::ZZ &x, NTL::ZZ &distance, NTL::ZZ &steps, long &index, bool &wild);

    inline NTL::ZZ setAlpha();
};

struct DistinguishedPoint {
    NTL::ZZ distance;
    bool isWild;
    long id;
};

#endif //LOGSPICKING_WEINERKANGAROO_H
