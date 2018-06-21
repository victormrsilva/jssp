NAME 
ROWS
 N  OBJ
 E  inicio(m1,t1)
 E  inicio(m2,t1)
 E  inicio(m3,t1)
 E  c29(1,2)
 E  c29(1,3)
 E  c29(1,4)
 E  c29(1,5)
 E  c29(1,6)
 E  c29(1,7)
 E  c29(1,8)
 E  c29(1,9)
 E  c29(1,10)
 E  c29(2,2)
 E  c29(2,3)
 E  c29(2,4)
 E  c29(2,5)
 E  c29(2,6)
 E  c29(2,7)
 E  c29(2,8)
 E  c29(2,9)
 E  c29(2,10)
 E  c29(3,2)
 E  c29(3,3)
 E  c29(3,4)
 E  c29(3,5)
 E  c29(3,6)
 E  c29(3,7)
 E  c29(3,8)
 E  c29(3,9)
 E  c29(3,10)
 E  espera_ini(m2,j1,t1)
 E  espera_ini(m1,j2,t1)
 E  espera_ini(m2,j3,t1)
 E  execute_wait(3,1,2)
 E  execute_wait(3,1,3)
 E  execute_wait(3,1,4)
 E  execute_wait(3,1,5)
 E  execute_wait(2,2,2)
 E  execute_wait(2,2,3)
 E  execute_wait(2,2,4)
 E  execute_wait(2,2,5)
 E  execute_wait(3,3,2)
 E  execute_wait(3,3,3)
 E  execute_wait(3,3,4)
 E  execute_wait(3,3,5)
 E  execute_wait(3,3,6)
 E  execute_wait(1,1,3)
 E  execute_wait(1,1,4)
 E  execute_wait(1,1,5)
 E  execute_wait(1,1,6)
 E  execute_wait(1,1,7)
 E  execute_wait(3,2,2)
 E  execute_wait(3,2,3)
 E  execute_wait(3,2,4)
 E  execute_wait(3,2,5)
 E  execute_wait(3,2,6)
 E  execute_wait(2,3,2)
 E  execute_wait(2,3,3)
 E  execute_wait(2,3,4)
 E  execute_wait(2,3,5)
 E  execute_wait(2,3,6)
 E  execute_wait(2,3,7)
 E  execute_wait(2,1,4)
 E  execute_wait(2,1,5)
 E  execute_wait(2,1,6)
 E  execute_wait(2,1,7)
 E  execute_wait(2,1,8)
 E  execute_wait(1,2,4)
 E  execute_wait(1,2,5)
 E  execute_wait(1,2,6)
 E  execute_wait(1,2,7)
 E  execute_wait(1,2,8)
 E  execute_wait(1,3,4)
 E  execute_wait(1,3,5)
 E  execute_wait(1,3,6)
 E  execute_wait(1,3,7)
 E  execute_wait(1,3,8)
 E  execute_wait(1,3,9)
 E  final_execute(3,1)
 E  final_execute(2,2)
 E  final_execute(3,3)
 E  final_execute(1,1)
 E  final_execute(3,2)
 E  final_execute(2,3)
 E  final_execute(2,1)
 E  final_execute(1,2)
 E  final_execute(1,3)
 G  makespan(1)
 G  makespan(2)
 G  makespan(3)
COLUMNS
    MARKER    'MARKER'                 'INTORG'
    f(1,1)    inicio(m1,t1)  -1
    f(1,1)    c29(1,2)  1
    f(1,2)    c29(1,2)  -1
    f(1,2)    c29(1,3)  1
    f(1,3)    c29(1,3)  -1
    f(1,3)    c29(1,4)  1
    x(1,1,3)  c29(1,3)  -1
    x(1,1,3)  c29(1,4)  1
    x(1,1,3)  execute_wait(1,1,3)  -1
    x(1,1,3)  execute_wait(2,1,4)  1
    e(1,1,3)  execute_wait(1,1,3)  -1
    e(1,1,3)  execute_wait(1,1,4)  1
    f(1,4)    c29(1,4)  -1
    f(1,4)    c29(1,5)  1
    x(1,1,4)  c29(1,4)  -1
    x(1,1,4)  c29(1,5)  1
    x(1,1,4)  execute_wait(1,1,4)  -1
    x(1,1,4)  execute_wait(2,1,5)  1
    e(1,1,4)  execute_wait(1,1,4)  -1
    e(1,1,4)  execute_wait(1,1,5)  1
    x(1,2,4)  c29(1,4)  -1
    x(1,2,4)  c29(1,6)  1
    x(1,2,4)  execute_wait(1,2,4)  -1
    x(1,2,4)  makespan(2)  -5
    e(1,2,4)  execute_wait(1,2,4)  -1
    e(1,2,4)  execute_wait(1,2,5)  1
    x(1,3,4)  c29(1,4)  -1
    x(1,3,4)  c29(1,5)  1
    x(1,3,4)  execute_wait(1,3,4)  -1
    x(1,3,4)  makespan(3)  -4
    e(1,3,4)  execute_wait(1,3,4)  -1
    e(1,3,4)  execute_wait(1,3,5)  1
    f(1,5)    c29(1,5)  -1
    f(1,5)    c29(1,6)  1
    x(1,1,5)  c29(1,5)  -1
    x(1,1,5)  c29(1,6)  1
    x(1,1,5)  execute_wait(1,1,5)  -1
    x(1,1,5)  execute_wait(2,1,6)  1
    e(1,1,5)  execute_wait(1,1,5)  -1
    e(1,1,5)  execute_wait(1,1,6)  1
    x(1,2,5)  c29(1,5)  -1
    x(1,2,5)  c29(1,7)  1
    x(1,2,5)  execute_wait(1,2,5)  -1
    x(1,2,5)  makespan(2)  -6
    e(1,2,5)  execute_wait(1,2,5)  -1
    e(1,2,5)  execute_wait(1,2,6)  1
    x(1,3,5)  c29(1,5)  -1
    x(1,3,5)  c29(1,6)  1
    x(1,3,5)  execute_wait(1,3,5)  -1
    x(1,3,5)  makespan(3)  -5
    e(1,3,5)  execute_wait(1,3,5)  -1
    e(1,3,5)  execute_wait(1,3,6)  1
    f(1,6)    c29(1,6)  -1
    f(1,6)    c29(1,7)  1
    x(1,1,6)  c29(1,6)  -1
    x(1,1,6)  c29(1,7)  1
    x(1,1,6)  execute_wait(1,1,6)  -1
    x(1,1,6)  execute_wait(2,1,7)  1
    e(1,1,6)  execute_wait(1,1,6)  -1
    e(1,1,6)  execute_wait(1,1,7)  1
    x(1,2,6)  c29(1,6)  -1
    x(1,2,6)  c29(1,8)  1
    x(1,2,6)  execute_wait(1,2,6)  -1
    x(1,2,6)  makespan(2)  -7
    e(1,2,6)  execute_wait(1,2,6)  -1
    e(1,2,6)  execute_wait(1,2,7)  1
    x(1,3,6)  c29(1,6)  -1
    x(1,3,6)  c29(1,7)  1
    x(1,3,6)  execute_wait(1,3,6)  -1
    x(1,3,6)  makespan(3)  -6
    e(1,3,6)  execute_wait(1,3,6)  -1
    e(1,3,6)  execute_wait(1,3,7)  1
    f(1,7)    c29(1,7)  -1
    f(1,7)    c29(1,8)  1
    x(1,1,7)  c29(1,7)  -1
    x(1,1,7)  c29(1,8)  1
    x(1,1,7)  execute_wait(1,1,7)  -1
    x(1,1,7)  execute_wait(2,1,8)  1
    e(1,1,7)  execute_wait(1,1,7)  -1
    e(1,1,7)  final_execute(1,1)  1
    x(1,2,7)  c29(1,7)  -1
    x(1,2,7)  c29(1,9)  1
    x(1,2,7)  execute_wait(1,2,7)  -1
    x(1,2,7)  makespan(2)  -8
    e(1,2,7)  execute_wait(1,2,7)  -1
    e(1,2,7)  execute_wait(1,2,8)  1
    x(1,3,7)  c29(1,7)  -1
    x(1,3,7)  c29(1,8)  1
    x(1,3,7)  execute_wait(1,3,7)  -1
    x(1,3,7)  makespan(3)  -7
    e(1,3,7)  execute_wait(1,3,7)  -1
    e(1,3,7)  execute_wait(1,3,8)  1
    f(1,8)    c29(1,8)  -1
    f(1,8)    c29(1,9)  1
    x(1,1,8)  c29(1,8)  -1
    x(1,1,8)  c29(1,9)  1
    x(1,1,8)  final_execute(1,1)  -1
    x(1,1,8)  final_execute(2,1)  1
    e(1,1,8)  OBJ       0
    x(1,2,8)  c29(1,8)  -1
    x(1,2,8)  c29(1,10)  1
    x(1,2,8)  execute_wait(1,2,8)  -1
    x(1,2,8)  makespan(2)  -9
    e(1,2,8)  execute_wait(1,2,8)  -1
    e(1,2,8)  final_execute(1,2)  1
    x(1,3,8)  c29(1,8)  -1
    x(1,3,8)  c29(1,9)  1
    x(1,3,8)  execute_wait(1,3,8)  -1
    x(1,3,8)  makespan(3)  -8
    e(1,3,8)  execute_wait(1,3,8)  -1
    e(1,3,8)  execute_wait(1,3,9)  1
    f(1,9)    c29(1,9)  -1
    f(1,9)    c29(1,10)  1
    x(1,2,9)  c29(1,9)  -1
    x(1,2,9)  final_execute(1,2)  -1
    x(1,2,9)  makespan(2)  -10
    e(1,2,9)  OBJ       0
    x(1,3,9)  c29(1,9)  -1
    x(1,3,9)  c29(1,10)  1
    x(1,3,9)  execute_wait(1,3,9)  -1
    x(1,3,9)  makespan(3)  -9
    e(1,3,9)  execute_wait(1,3,9)  -1
    e(1,3,9)  final_execute(1,3)  1
    f(1,10)   c29(1,10)  -1
    x(1,3,10)  c29(1,10)  -1
    x(1,3,10)  final_execute(1,3)  -1
    x(1,3,10)  makespan(3)  -10
    e(1,3,10)  OBJ       0
    f(2,1)    inicio(m2,t1)  -1
    f(2,1)    c29(2,2)  1
    x(2,2,1)  inicio(m2,t1)  -1
    x(2,2,1)  c29(2,2)  1
    x(2,2,1)  espera_ini(m1,j2,t1)  -1
    x(2,2,1)  execute_wait(3,2,2)  1
    e(2,2,1)  espera_ini(m1,j2,t1)  -1
    e(2,2,1)  execute_wait(2,2,2)  1
    f(2,2)    c29(2,2)  -1
    f(2,2)    c29(2,3)  1
    x(2,2,2)  c29(2,2)  -1
    x(2,2,2)  c29(2,3)  1
    x(2,2,2)  execute_wait(2,2,2)  -1
    x(2,2,2)  execute_wait(3,2,3)  1
    e(2,2,2)  execute_wait(2,2,2)  -1
    e(2,2,2)  execute_wait(2,2,3)  1
    x(2,3,2)  c29(2,2)  -1
    x(2,3,2)  c29(2,4)  1
    x(2,3,2)  execute_wait(2,3,2)  -1
    x(2,3,2)  execute_wait(1,3,4)  1
    e(2,3,2)  execute_wait(2,3,2)  -1
    e(2,3,2)  execute_wait(2,3,3)  1
    f(2,3)    c29(2,3)  -1
    f(2,3)    c29(2,4)  1
    x(2,2,3)  c29(2,3)  -1
    x(2,2,3)  c29(2,4)  1
    x(2,2,3)  execute_wait(2,2,3)  -1
    x(2,2,3)  execute_wait(3,2,4)  1
    e(2,2,3)  execute_wait(2,2,3)  -1
    e(2,2,3)  execute_wait(2,2,4)  1
    x(2,3,3)  c29(2,3)  -1
    x(2,3,3)  c29(2,5)  1
    x(2,3,3)  execute_wait(2,3,3)  -1
    x(2,3,3)  execute_wait(1,3,5)  1
    e(2,3,3)  execute_wait(2,3,3)  -1
    e(2,3,3)  execute_wait(2,3,4)  1
    f(2,4)    c29(2,4)  -1
    f(2,4)    c29(2,5)  1
    x(2,1,4)  c29(2,4)  -1
    x(2,1,4)  c29(2,6)  1
    x(2,1,4)  execute_wait(2,1,4)  -1
    x(2,1,4)  makespan(1)  -5
    e(2,1,4)  execute_wait(2,1,4)  -1
    e(2,1,4)  execute_wait(2,1,5)  1
    x(2,2,4)  c29(2,4)  -1
    x(2,2,4)  c29(2,5)  1
    x(2,2,4)  execute_wait(2,2,4)  -1
    x(2,2,4)  execute_wait(3,2,5)  1
    e(2,2,4)  execute_wait(2,2,4)  -1
    e(2,2,4)  execute_wait(2,2,5)  1
    x(2,3,4)  c29(2,4)  -1
    x(2,3,4)  c29(2,6)  1
    x(2,3,4)  execute_wait(2,3,4)  -1
    x(2,3,4)  execute_wait(1,3,6)  1
    e(2,3,4)  execute_wait(2,3,4)  -1
    e(2,3,4)  execute_wait(2,3,5)  1
    f(2,5)    c29(2,5)  -1
    f(2,5)    c29(2,6)  1
    x(2,1,5)  c29(2,5)  -1
    x(2,1,5)  c29(2,7)  1
    x(2,1,5)  execute_wait(2,1,5)  -1
    x(2,1,5)  makespan(1)  -6
    e(2,1,5)  execute_wait(2,1,5)  -1
    e(2,1,5)  execute_wait(2,1,6)  1
    x(2,2,5)  c29(2,5)  -1
    x(2,2,5)  c29(2,6)  1
    x(2,2,5)  execute_wait(2,2,5)  -1
    x(2,2,5)  execute_wait(3,2,6)  1
    e(2,2,5)  execute_wait(2,2,5)  -1
    e(2,2,5)  final_execute(2,2)  1
    x(2,3,5)  c29(2,5)  -1
    x(2,3,5)  c29(2,7)  1
    x(2,3,5)  execute_wait(2,3,5)  -1
    x(2,3,5)  execute_wait(1,3,7)  1
    e(2,3,5)  execute_wait(2,3,5)  -1
    e(2,3,5)  execute_wait(2,3,6)  1
    f(2,6)    c29(2,6)  -1
    f(2,6)    c29(2,7)  1
    x(2,1,6)  c29(2,6)  -1
    x(2,1,6)  c29(2,8)  1
    x(2,1,6)  execute_wait(2,1,6)  -1
    x(2,1,6)  makespan(1)  -7
    e(2,1,6)  execute_wait(2,1,6)  -1
    e(2,1,6)  execute_wait(2,1,7)  1
    x(2,2,6)  c29(2,6)  -1
    x(2,2,6)  c29(2,7)  1
    x(2,2,6)  final_execute(2,2)  -1
    x(2,2,6)  final_execute(3,2)  1
    e(2,2,6)  OBJ       0
    x(2,3,6)  c29(2,6)  -1
    x(2,3,6)  c29(2,8)  1
    x(2,3,6)  execute_wait(2,3,6)  -1
    x(2,3,6)  execute_wait(1,3,8)  1
    e(2,3,6)  execute_wait(2,3,6)  -1
    e(2,3,6)  execute_wait(2,3,7)  1
    f(2,7)    c29(2,7)  -1
    f(2,7)    c29(2,8)  1
    x(2,1,7)  c29(2,7)  -1
    x(2,1,7)  c29(2,9)  1
    x(2,1,7)  execute_wait(2,1,7)  -1
    x(2,1,7)  makespan(1)  -8
    e(2,1,7)  execute_wait(2,1,7)  -1
    e(2,1,7)  execute_wait(2,1,8)  1
    x(2,3,7)  c29(2,7)  -1
    x(2,3,7)  c29(2,9)  1
    x(2,3,7)  execute_wait(2,3,7)  -1
    x(2,3,7)  execute_wait(1,3,9)  1
    e(2,3,7)  execute_wait(2,3,7)  -1
    e(2,3,7)  final_execute(2,3)  1
    f(2,8)    c29(2,8)  -1
    f(2,8)    c29(2,9)  1
    x(2,1,8)  c29(2,8)  -1
    x(2,1,8)  c29(2,10)  1
    x(2,1,8)  execute_wait(2,1,8)  -1
    x(2,1,8)  makespan(1)  -9
    e(2,1,8)  execute_wait(2,1,8)  -1
    e(2,1,8)  final_execute(2,1)  1
    x(2,3,8)  c29(2,8)  -1
    x(2,3,8)  c29(2,10)  1
    x(2,3,8)  final_execute(2,3)  -1
    x(2,3,8)  final_execute(1,3)  1
    e(2,3,8)  OBJ       0
    f(2,9)    c29(2,9)  -1
    f(2,9)    c29(2,10)  1
    x(2,1,9)  c29(2,9)  -1
    x(2,1,9)  final_execute(2,1)  -1
    x(2,1,9)  makespan(1)  -10
    e(2,1,9)  OBJ       0
    f(2,10)   c29(2,10)  -1
    f(3,1)    inicio(m3,t1)  -1
    f(3,1)    c29(3,2)  1
    x(3,1,1)  inicio(m3,t1)  -1
    x(3,1,1)  c29(3,3)  1
    x(3,1,1)  espera_ini(m2,j1,t1)  -1
    x(3,1,1)  execute_wait(1,1,3)  1
    e(3,1,1)  espera_ini(m2,j1,t1)  -1
    e(3,1,1)  execute_wait(3,1,2)  1
    x(3,3,1)  inicio(m3,t1)  -1
    x(3,3,1)  c29(3,2)  1
    x(3,3,1)  espera_ini(m2,j3,t1)  -1
    x(3,3,1)  execute_wait(2,3,2)  1
    e(3,3,1)  espera_ini(m2,j3,t1)  -1
    e(3,3,1)  execute_wait(3,3,2)  1
    f(3,2)    c29(3,2)  -1
    f(3,2)    c29(3,3)  1
    x(3,1,2)  c29(3,2)  -1
    x(3,1,2)  c29(3,4)  1
    x(3,1,2)  execute_wait(3,1,2)  -1
    x(3,1,2)  execute_wait(1,1,4)  1
    e(3,1,2)  execute_wait(3,1,2)  -1
    e(3,1,2)  execute_wait(3,1,3)  1
    x(3,2,2)  c29(3,2)  -1
    x(3,2,2)  c29(3,4)  1
    x(3,2,2)  execute_wait(3,2,2)  -1
    x(3,2,2)  execute_wait(1,2,4)  1
    e(3,2,2)  execute_wait(3,2,2)  -1
    e(3,2,2)  execute_wait(3,2,3)  1
    x(3,3,2)  c29(3,2)  -1
    x(3,3,2)  c29(3,3)  1
    x(3,3,2)  execute_wait(3,3,2)  -1
    x(3,3,2)  execute_wait(2,3,3)  1
    e(3,3,2)  execute_wait(3,3,2)  -1
    e(3,3,2)  execute_wait(3,3,3)  1
    f(3,3)    c29(3,3)  -1
    f(3,3)    c29(3,4)  1
    x(3,1,3)  c29(3,3)  -1
    x(3,1,3)  c29(3,5)  1
    x(3,1,3)  execute_wait(3,1,3)  -1
    x(3,1,3)  execute_wait(1,1,5)  1
    e(3,1,3)  execute_wait(3,1,3)  -1
    e(3,1,3)  execute_wait(3,1,4)  1
    x(3,2,3)  c29(3,3)  -1
    x(3,2,3)  c29(3,5)  1
    x(3,2,3)  execute_wait(3,2,3)  -1
    x(3,2,3)  execute_wait(1,2,5)  1
    e(3,2,3)  execute_wait(3,2,3)  -1
    e(3,2,3)  execute_wait(3,2,4)  1
    x(3,3,3)  c29(3,3)  -1
    x(3,3,3)  c29(3,4)  1
    x(3,3,3)  execute_wait(3,3,3)  -1
    x(3,3,3)  execute_wait(2,3,4)  1
    e(3,3,3)  execute_wait(3,3,3)  -1
    e(3,3,3)  execute_wait(3,3,4)  1
    f(3,4)    c29(3,4)  -1
    f(3,4)    c29(3,5)  1
    x(3,1,4)  c29(3,4)  -1
    x(3,1,4)  c29(3,6)  1
    x(3,1,4)  execute_wait(3,1,4)  -1
    x(3,1,4)  execute_wait(1,1,6)  1
    e(3,1,4)  execute_wait(3,1,4)  -1
    e(3,1,4)  execute_wait(3,1,5)  1
    x(3,2,4)  c29(3,4)  -1
    x(3,2,4)  c29(3,6)  1
    x(3,2,4)  execute_wait(3,2,4)  -1
    x(3,2,4)  execute_wait(1,2,6)  1
    e(3,2,4)  execute_wait(3,2,4)  -1
    e(3,2,4)  execute_wait(3,2,5)  1
    x(3,3,4)  c29(3,4)  -1
    x(3,3,4)  c29(3,5)  1
    x(3,3,4)  execute_wait(3,3,4)  -1
    x(3,3,4)  execute_wait(2,3,5)  1
    e(3,3,4)  execute_wait(3,3,4)  -1
    e(3,3,4)  execute_wait(3,3,5)  1
    f(3,5)    c29(3,5)  -1
    f(3,5)    c29(3,6)  1
    x(3,1,5)  c29(3,5)  -1
    x(3,1,5)  c29(3,7)  1
    x(3,1,5)  execute_wait(3,1,5)  -1
    x(3,1,5)  execute_wait(1,1,7)  1
    e(3,1,5)  execute_wait(3,1,5)  -1
    e(3,1,5)  final_execute(3,1)  1
    x(3,2,5)  c29(3,5)  -1
    x(3,2,5)  c29(3,7)  1
    x(3,2,5)  execute_wait(3,2,5)  -1
    x(3,2,5)  execute_wait(1,2,7)  1
    e(3,2,5)  execute_wait(3,2,5)  -1
    e(3,2,5)  execute_wait(3,2,6)  1
    x(3,3,5)  c29(3,5)  -1
    x(3,3,5)  c29(3,6)  1
    x(3,3,5)  execute_wait(3,3,5)  -1
    x(3,3,5)  execute_wait(2,3,6)  1
    e(3,3,5)  execute_wait(3,3,5)  -1
    e(3,3,5)  execute_wait(3,3,6)  1
    f(3,6)    c29(3,6)  -1
    f(3,6)    c29(3,7)  1
    x(3,1,6)  c29(3,6)  -1
    x(3,1,6)  c29(3,8)  1
    x(3,1,6)  final_execute(3,1)  -1
    x(3,1,6)  final_execute(1,1)  1
    e(3,1,6)  OBJ       0
    x(3,2,6)  c29(3,6)  -1
    x(3,2,6)  c29(3,8)  1
    x(3,2,6)  execute_wait(3,2,6)  -1
    x(3,2,6)  execute_wait(1,2,8)  1
    e(3,2,6)  execute_wait(3,2,6)  -1
    e(3,2,6)  final_execute(3,2)  1
    x(3,3,6)  c29(3,6)  -1
    x(3,3,6)  c29(3,7)  1
    x(3,3,6)  execute_wait(3,3,6)  -1
    x(3,3,6)  execute_wait(2,3,7)  1
    e(3,3,6)  execute_wait(3,3,6)  -1
    e(3,3,6)  final_execute(3,3)  1
    f(3,7)    c29(3,7)  -1
    f(3,7)    c29(3,8)  1
    x(3,2,7)  c29(3,7)  -1
    x(3,2,7)  c29(3,9)  1
    x(3,2,7)  final_execute(3,2)  -1
    x(3,2,7)  final_execute(1,2)  1
    e(3,2,7)  OBJ       0
    x(3,3,7)  c29(3,7)  -1
    x(3,3,7)  c29(3,8)  1
    x(3,3,7)  final_execute(3,3)  -1
    x(3,3,7)  final_execute(2,3)  1
    e(3,3,7)  OBJ       0
    f(3,8)    c29(3,8)  -1
    f(3,8)    c29(3,9)  1
    f(3,9)    c29(3,9)  -1
    f(3,9)    c29(3,10)  1
    f(3,10)   c29(3,10)  -1
    C         OBJ       1
    C         makespan(1)  1
    C         makespan(2)  1
    C         makespan(3)  1
    MARKER    'MARKER'                 'INTEND'
RHS
    RHS1      inicio(m1,t1)  -1
    RHS1      inicio(m2,t1)  -1
    RHS1      inicio(m3,t1)  -1
    RHS1      espera_ini(m2,j1,t1)  -1
    RHS1      espera_ini(m1,j2,t1)  -1
    RHS1      espera_ini(m2,j3,t1)  -1
BOUNDS
 BV BND1      f(1,1)  
 BV BND1      f(1,2)  
 BV BND1      f(1,3)  
 BV BND1      x(1,1,3)
 BV BND1      e(1,1,3)
 BV BND1      f(1,4)  
 BV BND1      x(1,1,4)
 BV BND1      e(1,1,4)
 BV BND1      x(1,2,4)
 BV BND1      e(1,2,4)
 BV BND1      x(1,3,4)
 BV BND1      e(1,3,4)
 BV BND1      f(1,5)  
 BV BND1      x(1,1,5)
 BV BND1      e(1,1,5)
 BV BND1      x(1,2,5)
 BV BND1      e(1,2,5)
 BV BND1      x(1,3,5)
 BV BND1      e(1,3,5)
 BV BND1      f(1,6)  
 BV BND1      x(1,1,6)
 BV BND1      e(1,1,6)
 BV BND1      x(1,2,6)
 BV BND1      e(1,2,6)
 BV BND1      x(1,3,6)
 BV BND1      e(1,3,6)
 BV BND1      f(1,7)  
 BV BND1      x(1,1,7)
 BV BND1      e(1,1,7)
 BV BND1      x(1,2,7)
 BV BND1      e(1,2,7)
 BV BND1      x(1,3,7)
 BV BND1      e(1,3,7)
 BV BND1      f(1,8)  
 BV BND1      x(1,1,8)
 BV BND1      e(1,1,8)
 BV BND1      x(1,2,8)
 BV BND1      e(1,2,8)
 BV BND1      x(1,3,8)
 BV BND1      e(1,3,8)
 BV BND1      f(1,9)  
 BV BND1      x(1,2,9)
 BV BND1      e(1,2,9)
 BV BND1      x(1,3,9)
 BV BND1      e(1,3,9)
 BV BND1      f(1,10) 
 BV BND1      x(1,3,10)
 BV BND1      e(1,3,10)
 BV BND1      f(2,1)  
 BV BND1      x(2,2,1)
 BV BND1      e(2,2,1)
 BV BND1      f(2,2)  
 BV BND1      x(2,2,2)
 BV BND1      e(2,2,2)
 BV BND1      x(2,3,2)
 BV BND1      e(2,3,2)
 BV BND1      f(2,3)  
 BV BND1      x(2,2,3)
 BV BND1      e(2,2,3)
 BV BND1      x(2,3,3)
 BV BND1      e(2,3,3)
 BV BND1      f(2,4)  
 BV BND1      x(2,1,4)
 BV BND1      e(2,1,4)
 BV BND1      x(2,2,4)
 BV BND1      e(2,2,4)
 BV BND1      x(2,3,4)
 BV BND1      e(2,3,4)
 BV BND1      f(2,5)  
 BV BND1      x(2,1,5)
 BV BND1      e(2,1,5)
 BV BND1      x(2,2,5)
 BV BND1      e(2,2,5)
 BV BND1      x(2,3,5)
 BV BND1      e(2,3,5)
 BV BND1      f(2,6)  
 BV BND1      x(2,1,6)
 BV BND1      e(2,1,6)
 BV BND1      x(2,2,6)
 BV BND1      e(2,2,6)
 BV BND1      x(2,3,6)
 BV BND1      e(2,3,6)
 BV BND1      f(2,7)  
 BV BND1      x(2,1,7)
 BV BND1      e(2,1,7)
 BV BND1      x(2,3,7)
 BV BND1      e(2,3,7)
 BV BND1      f(2,8)  
 BV BND1      x(2,1,8)
 BV BND1      e(2,1,8)
 BV BND1      x(2,3,8)
 BV BND1      e(2,3,8)
 BV BND1      f(2,9)  
 BV BND1      x(2,1,9)
 BV BND1      e(2,1,9)
 BV BND1      f(2,10) 
 BV BND1      f(3,1)  
 BV BND1      x(3,1,1)
 BV BND1      e(3,1,1)
 BV BND1      x(3,3,1)
 BV BND1      e(3,3,1)
 BV BND1      f(3,2)  
 BV BND1      x(3,1,2)
 BV BND1      e(3,1,2)
 BV BND1      x(3,2,2)
 BV BND1      e(3,2,2)
 BV BND1      x(3,3,2)
 BV BND1      e(3,3,2)
 BV BND1      f(3,3)  
 BV BND1      x(3,1,3)
 BV BND1      e(3,1,3)
 BV BND1      x(3,2,3)
 BV BND1      e(3,2,3)
 BV BND1      x(3,3,3)
 BV BND1      e(3,3,3)
 BV BND1      f(3,4)  
 BV BND1      x(3,1,4)
 BV BND1      e(3,1,4)
 BV BND1      x(3,2,4)
 BV BND1      e(3,2,4)
 BV BND1      x(3,3,4)
 BV BND1      e(3,3,4)
 BV BND1      f(3,5)  
 BV BND1      x(3,1,5)
 BV BND1      e(3,1,5)
 BV BND1      x(3,2,5)
 BV BND1      e(3,2,5)
 BV BND1      x(3,3,5)
 BV BND1      e(3,3,5)
 BV BND1      f(3,6)  
 BV BND1      x(3,1,6)
 BV BND1      e(3,1,6)
 BV BND1      x(3,2,6)
 BV BND1      e(3,2,6)
 BV BND1      x(3,3,6)
 BV BND1      e(3,3,6)
 BV BND1      f(3,7)  
 BV BND1      x(3,2,7)
 BV BND1      e(3,2,7)
 BV BND1      x(3,3,7)
 BV BND1      e(3,3,7)
 BV BND1      f(3,8)  
 BV BND1      f(3,9)  
 BV BND1      f(3,10) 
 LI BND1      C         0
ENDATA
