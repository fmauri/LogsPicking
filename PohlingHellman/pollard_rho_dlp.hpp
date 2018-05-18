#ifndef POLLARD_RHO_DLP_INCLUDED
#define POLLARD_RHO_DLP_INCLUDED

#include "../Support/utilities.hpp"
#include <NTL/ZZ.h>
#include <omp.h>
#include <algorithm>
#include <unordered_map>

const long THETA = 65536; // 2^16

struct DistinguishedPoint {
    NTL::ZZ Position;
    NTL::ZZ ExponentA;
    NTL::ZZ ExponentB;
};

NTL::ZZ PollardRhoDLP(NTL::ZZ p, NTL::ZZ q, NTL::ZZ a, NTL::ZZ b);

#endif