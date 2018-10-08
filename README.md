# jssp
Job-shop sheduling problem

compilation: cmake CMakeLists.txt & make
OBS: Needs Gurobi 8.0 installed in /opt/gurobi800/linux64. If it's installed elsewhere, change set(LP_COMPILE_FLAGS "-DGRB -I/opt/gurobi800/linux64/include") line in CMakeLists.txt

Execution:
jssp instanceFile timeLimit execute formulation
timeLimit: a integer with limited time of execution. -1 will be maximum time needed
execute: 1 for execute the lp. 0 only generates lp
formulation: 
F for Flow (packing formulation)
C for Compact (BigM) 
K for Kondilli 
Fe for Fernando (machine formulation) 
