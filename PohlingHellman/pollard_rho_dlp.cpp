#include "pollard_rho_dlp.hpp"

inline int Partition(NTL::ZZ x)
{
    return x % 3;
}

inline bool IsDistinguishedPoint(NTL::ZZ x)
{
    return x % THETA == 0;
}

void PerformStep(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b, NTL::ZZ &step, NTL::ZZ g, NTL::ZZ y, NTL::ZZ p, NTL::ZZ q)
{
    switch (Partition(x))
    {
        case 0:
            a = NTL::MulMod(a, 2, q);
            b = NTL::MulMod(b, 2, q);
            x = NTL::MulMod(x, x, p);
            break;
        case 1:
            a = NTL::AddMod(a, 1, q);
            x = NTL::MulMod(x, g, p);
            break;
        case 2:
            b = NTL::AddMod(b, 1, q);
            x = NTL::MulMod(x, y, p);
            break;
        default:
            throw std::invalid_argument("Something went horribly wrong.");
    }
 
    step += 1;
}

void InitializeStartingPoint(NTL::ZZ &x, NTL::ZZ &a, NTL::ZZ &b, NTL::ZZ &step, NTL::ZZ g, NTL::ZZ y, NTL::ZZ p)
{
    a = FindRandomElementOfGroup(p);
    b = FindRandomElementOfGroup(p);
    x = NTL::MulMod(NTL::PowerMod(g, a, p), NTL::PowerMod(y, b, p), p);
    step = NTL::ZZ(0);
}

NTL::ZZ PollardRhoDLP(NTL::ZZ p, NTL::ZZ q, NTL::ZZ a, NTL::ZZ b)
{
    std::unordered_map<std::string, DistinguishedPoint> points;
    NTL::ZZ pX, pA, pB, pStep, res;

    #pragma omp parallel for shared(points, res) private(pX, pA, pB, pStep)
    for (int i = 0; i < omp_get_num_threads(); ++i)
    {
        while (res == 0)
        {
            InitializeStartingPoint(pX, pA, pB, pStep, a, b, p);

            while (!IsDistinguishedPoint(pX))
            {
                if (res != 0 || pStep > THETA * 20) break;

                PerformStep(pX, pA, pB, pStep, a, b, p, q);
            }

            if (res != 0 || pStep > THETA * 20) continue;

            std::string id = ZzToString(pX);
            auto collision = points.find(id);

            if (collision != points.end())
            {
                NTL::ZZ pC = collision->second.ExponentA;
                NTL::ZZ pD = collision->second.ExponentB;
                
                NTL::ZZ r = NTL::SubMod(pD, pB, q);
                if (r == 0) continue;
                
                NTL::ZZ candidate = NTL::MulMod(NTL::InvMod(r, q), NTL::SubMod(pA, pC, q), q);
                if (NTL::PowerMod(a, candidate, p) != b) continue;
                
                res = candidate;
            }
            else
            {
                points[id] = { pX, pA, pB };
            }
        }
    }
    
    return res;
}