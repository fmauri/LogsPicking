#include <iostream>
#include <omp.h>
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
//    if (pollarRho.isBeta(result)) {
//        std::cout << "Correct\n";
//    } else {
//        std::cout << "Wrong\n";
//    }
//    /*
//     * Kangaroo Weiner
//     */
//    std::cout << "\nKangaroo Weiner\n";
//    WeinerKangaroo weinerKangaroo;
//    weinerKangaroo.printLog();
//    result = weinerKangaroo.searchCollisions();
//    std::cout << "x = " << result << std::endl;
//    if (weinerKangaroo.isBeta(result)) {
//        std::cout << "Correct\n";
//    } else {
//        std::cout << "Wrong\n";
//    }
    /*
     * PohlingHellman
     */
    PohlingHellman pohlingHellman;
    omp_set_num_threads(8);
    std::cout << "Computing challenge" << std::endl;
    pohlingHellman.GenerateChallenge();
    std::cout << "Done\n" << "Goal is:\n" << pohlingHellman.getGoal() << "\n";
    std::cout << "Computing Pohling-Hellman" << std::endl;
    result = pohlingHellman.calculateReduction();
    std::cout << "Calculated\n" << result << std::endl;

    return 0;
}