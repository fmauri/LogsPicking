#include "utilities.hpp"

NTL::ZZ FindStrongPrime(int l)
{
    NTL::ZZ q = NTL::GenGermainPrime_ZZ(l);
    return 2 * q + 1;
}

NTL::ZZ FindRandomElementOfGroup(NTL::ZZ p)
{
    return NTL::RandomBnd(p - 1) + NTL::ZZ(1);
}

NTL::ZZ FindGroupGenerator(NTL::ZZ p)
{
    while (true)
    {
        NTL::ZZ g = FindRandomElementOfGroup(p);
        g = NTL::MulMod(g, g, p);
        
        if (g == 1) continue;
        
        return g;
    }
}

std::string ZzToString(NTL::ZZ n)
{
    std::stringstream buffer;
    buffer << n;

    return buffer.str();
}