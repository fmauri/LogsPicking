#include <iostream>
#include "PollardRho/PollarRho.h"
#include "WeinerKangaroo/WeinerKangaroo.h"

int main() {
    NTL::ZZ result;
    /*
     * PollarRho
     */
    std::cout << "Pollard Rho\n";
    PollarRho pollarRho;
    pollarRho.printLog();
    result = pollarRho.searchXParallelPollard();
    std::cout << "x = " << result << std::endl;
    /*
     * Kangaroo Weiner
     */
    std::cout << "\nKangaroo Weiner\n";
    WeinerKangaroo weinerKangaroo;
    weinerKangaroo.printLog();
    result = weinerKangaroo.searchCollisions();
    std::cout << "x = " << result << std::endl;
    return 0;
}