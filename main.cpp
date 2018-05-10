#include <iostream>
#include "PollardRho/PollarRho.h"
#include "WeinerKangaroo/WeinerKangaroo.h"

int main() {
    NTL::ZZ result;
    /*
     * PollarRho
     */
    PollarRho pollarRho;
    pollarRho.printLog();
    result = pollarRho.searchXParallelPollard();
    std::cout << "x = " << result;
    /*
     * Kangaroo Pollar Rho
     */
    WeinerKangaroo weinerKangaroo;
    result = weinerKangaroo.searchCollisions();
    return 0;
}