cmake_minimum_required(VERSION 3.5)
project(CGraph)

set(CMAKE_CXX_STANDARD 14)

set(LP_COMPILE_FLAGS "-DCBC -I/opt/cbc-default/include/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LP_COMPILE_FLAGS} -fsanitize=address")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LP_COMPILE_FLAGS}")

set(LP_LINK_FLAGS "-L/opt/cbc-default/lib/ -lCbcSolver -lCbc -lCgl -lOsiClp -lClpSolver -lClp -lOsi -lCoinUtils -lcoinlapack -lcoinblas -lbz2 -lz -lgfortran -lgomp -lquadmath -ldl -lpthread")

add_executable(CGraph main.cpp lp.h lp.cpp build_cgraph.h build_cgraph.cpp cgraph.h cgraph.c clique.h clique.c macros.h
               node_heap.h node_heap.c vint_set.h vint_set.c vectormgm.h vectormgm.c strutils.h strutils.c memory.h memory.c)

target_link_libraries(CGraph ${LP_LINK_FLAGS})
