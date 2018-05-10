cmake_minimum_required(VERSION 3.10)
project(LogsPicking)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

add_executable(LogsPicking main.cpp PollardRho/PollarRho.cpp PollardRho/PollarRho.h DiscreteLog.h WeinerKangaroo/WeinerKangaroo.cpp WeinerKangaroo/WeinerKangaroo.h)

target_link_libraries(LogsPicking ntl m gmp)