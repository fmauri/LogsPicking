#include <iostream>
#include "PollardRho/PollardRho.h"
#include "WeinerKangaroo/WeinerKangaroo.h"
#include "PohlingHellman/PohlingHellman.h"

int main() {
    NTL::ZZ result;
    /*
     * PollardRho
     */
//    std::cout << "Pollard Rho\n";
//    PollardRho pollarRho;
//    pollarRho.printLog();
//    result = pollarRho.searchXParallelPollard();
//    std::cout << "x = " << result << std::endl;
//    /*
//     * Kangaroo Weiner
//     */
//    std::cout << "\nKangaroo Weiner\n";
//    WeinerKangaroo weinerKangaroo;
//    weinerKangaroo.printLog();
//    result = weinerKangaroo.searchCollisions();
//    std::cout << "x = " << result << std::endl;
    /*
     * PohlingHellman
     */
    std::cout << "\nPohling Hellman\n";
    PohlingHellman pohlingHellman;
    pohlingHellman.printLog();
    result = pohlingHellman.searchResult();
    std::cout << "x = " << result << std::endl;
    return 0;
}