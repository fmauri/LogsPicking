#ifndef UTILITIES_INCLUDED
#define UTILITIES_INCLUDED

#include <NTL/ZZ.h>
#include <sstream>
#include <string>

NTL::ZZ FindStrongPrime(int l);
NTL::ZZ FindRandomElementOfGroup(NTL::ZZ p);
NTL::ZZ FindGroupGenerator(NTL::ZZ p);
std::string ZzToString(NTL::ZZ n);

#endif