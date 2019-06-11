cmake_minimum_required(VERSION 3.5)
project(jssp)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 99)
#-Wall -g3  -std=c++11  
set(LP_COMPILE_FLAGS "-DGRB -I/opt/gurobi811/linux64/include")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -Wall -g3  ${LP_COMPILE_FLAGS}")
#-fsanitize=address
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LP_COMPILE_FLAGS}")

set(LP_LINK_FLAGS "-L/opt/gurobi811/linux64/lib -lgurobi81 -lgurobi_c++")
#set(LP_LINK_FLAGS "-L/opt/cbc-default/lib/ -lCbcSolver -lCbc -lCgl -lOsiClp -lClpSolver -lClp -lOsi -lCoinUtils -lcoinlapack -lcoinblas -lbz2 -lz -lgfortran -lgomp -lquadmath -ldl -lpthread")
#src/Flow.cpp src/Flow.hpp src/lp.cpp src/lp.h
add_executable(jssp main.cpp src/Gera.cpp src/Gera.hpp src/Compact.cpp src/Compact.hpp src/Fernando.cpp src/Fernando.hpp  src/Flow_testes.cpp src/Flow_testes.hpp 
               src/Instance.cpp src/Instance.hpp src/Kondili.cpp src/Kondili.hpp src/lp.hpp src/Callback.cpp src/Callback.hpp   src/Flow.cpp src/Flow.hpp
               src/cgraph/build_cgraph.cpp src/cgraph/build_cgraph.h src/cgraph/cgraph.c src/cgraph/cgraph.h src/cgraph/clique.c src/cgraph/clique.h 
               src/cgraph/lp.cpp src/cgraph/lp.h src/cgraph/macros.h src/cgraph/memory.c src/cgraph/memory.h src/cgraph/node_heap.c src/cgraph/node_heap.h 
               src/cgraph/BKGraph.cpp src/cgraph/BKGraph.hpp src/cgraph/BKVertex.hpp src/cgraph/BKVertex.cpp src/cgraph/bron_kerbosch.cpp src/cgraph/bron_kerbosch.h 
               src/cgraph/clique_extender.c src/cgraph/clique_extender.h src/cgraph/clique_merge.h src/cgraph/clique_merge.c src/cgraph/clique_separation.c src/cgraph/clique_separation.h 
               src/cgraph/strutils.c src/cgraph/strutils.h src/cgraph/vectormgm.c src/cgraph/vectormgm.h src/cgraph/vint_set.c src/cgraph/vint_set.h src/cgraph/vint_queue.c src/cgraph/vint_queue.h
               src/cgraph/spaths.h src/cgraph/spaths.c src/cgraph/oddhs.h src/cgraph/oddhs.c src/cgraph/digraph.h src/cgraph/cut.h src/cgraph/cut.c src/Hash.cpp)

target_link_libraries(jssp ${LP_LINK_FLAGS})
