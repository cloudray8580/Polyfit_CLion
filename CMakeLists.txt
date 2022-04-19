cmake_minimum_required(VERSION 3.16)
project(Polyfit_CLion)

set(CMAKE_CXX_STANDARD 17)

add_executable(Polyfit_CLion main.cpp utils.h Polyfit_1D.h MethodTesting.h S2Sampling.h FITingTree.h RMI.h EntropyHistogram.h Polyfit_2D.h SimpleRTree.h PGM.h experiments.h DBest.h VerdictDB.h Polyfit.h test.h RadixSpline.h Polyfit2D.h SWAB.h MultiResolutionTree.h ICTAI1997.h aRtree.h)

add_definitions(-DILOUSESTL -DIL_STD)
# ILOUSESTL: iostream.h: No such file or directory
# IL_STD: cerr was not declared in this scope; did you mean std::cerr

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

include_directories(
        /mnt/d/PGM-index/include
        /mnt/d/stxbtree/stx-btree-0.9/include
        /mnt/d/boostlinux/boost_1_74_0
        /mnt/d/CPLEXLinux/concert/include
        /mnt/d/CPLEXLinux/cplex/include
        /mnt/d/ensmallen/ensmallen-2.14.0/ensmallen-2.14.0/include
        /mnt/d/xgboost/include
        /mnt/d/RadixSpline/include
) # WSL file system format, if using cygwin: start from /cygdrive/d/

# remember to use the linux version library
find_library(ILOCPLEX ilocplex /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/)
find_library(CPLEX cplex /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/)
find_library(CONCERT concert /mnt/d/CPLEXLinux/concert/lib/x86-64_linux/static_pic/)
find_library(XGBOOST xgboost /mnt/d/xgboost/lib/)

target_link_libraries(Polyfit_CLion PUBLIC
        ${ILOCPLEX}
        ${CONCERT}
        ${CPLEX}
        libm.so # should use dynamic lib
        libpthread.so # should use dynamic lib
        libboost_serialization.so # link libboost_serialization.so at /usr/local/
        libmlpack.so # link libmlpack.so at /usr/lib/x86_64-linux-gnu/
        ${XGBOOST}
        ${CMAKE_DL_LIBS}
        )