NAME 
ROWS
 N  OBJ
 E  inicio_maquina(m1,t1)
 E  inicio_maquina(m2,t1)
 E  inicio_maquina(m3,t1)
 E  fluxo_maquina(1,2)
 E  fluxo_maquina(1,3)
 E  fluxo_maquina(1,4)
 E  fluxo_maquina(1,5)
 E  fluxo_maquina(1,6)
 E  fluxo_maquina(1,7)
 E  fluxo_maquina(1,8)
 E  fluxo_maquina(1,9)
 E  fluxo_maquina(1,10)
 E  fluxo_maquina(1,11)
 E  fluxo_maquina(1,12)
 E  fluxo_maquina(2,2)
 E  fluxo_maquina(2,3)
 E  fluxo_maquina(2,4)
 E  fluxo_maquina(2,5)
 E  fluxo_maquina(2,6)
 E  fluxo_maquina(2,7)
 E  fluxo_maquina(2,8)
 E  fluxo_maquina(2,9)
 E  fluxo_maquina(2,10)
 E  fluxo_maquina(2,11)
 E  fluxo_maquina(2,12)
 E  fluxo_maquina(3,2)
 E  fluxo_maquina(3,3)
 E  fluxo_maquina(3,4)
 E  fluxo_maquina(3,5)
 E  fluxo_maquina(3,6)
 E  fluxo_maquina(3,7)
 E  fluxo_maquina(3,8)
 E  fluxo_maquina(3,9)
 E  fluxo_maquina(3,10)
 E  fluxo_maquina(3,11)
 E  fluxo_maquina(3,12)
 E  inicio_espera(m0,j1,t1)
 E  inicio_espera(m0,j2,t1)
 E  inicio_espera(m1,j3,t1)
 E  fluxo_espera(1,1,2)
 E  fluxo_espera(1,1,3)
 E  fluxo_espera(1,1,4)
 E  fluxo_espera(1,1,5)
 E  fluxo_espera(1,1,6)
 E  fluxo_espera(1,1,7)
 E  fluxo_espera(1,2,2)
 E  fluxo_espera(1,2,3)
 E  fluxo_espera(1,2,4)
 E  fluxo_espera(1,2,5)
 E  fluxo_espera(1,2,6)
 E  fluxo_espera(1,2,7)
 E  fluxo_espera(2,3,2)
 E  fluxo_espera(2,3,3)
 E  fluxo_espera(2,3,4)
 E  fluxo_espera(2,3,5)
 E  fluxo_espera(2,1,2)
 E  fluxo_espera(2,1,3)
 E  fluxo_espera(2,1,4)
 E  fluxo_espera(2,1,5)
 E  fluxo_espera(2,1,6)
 E  fluxo_espera(2,1,7)
 E  fluxo_espera(2,1,8)
 E  fluxo_espera(3,2,2)
 E  fluxo_espera(3,2,3)
 E  fluxo_espera(3,2,4)
 E  fluxo_espera(3,2,5)
 E  fluxo_espera(3,2,6)
 E  fluxo_espera(3,2,7)
 E  fluxo_espera(3,2,8)
 E  fluxo_espera(1,3,5)
 E  fluxo_espera(1,3,6)
 E  fluxo_espera(1,3,7)
 E  fluxo_espera(1,3,8)
 E  fluxo_espera(1,3,9)
 E  fluxo_espera(3,1,4)
 E  fluxo_espera(3,1,5)
 E  fluxo_espera(3,1,6)
 E  fluxo_espera(3,1,7)
 E  fluxo_espera(3,1,8)
 E  fluxo_espera(3,1,9)
 E  fluxo_espera(3,1,10)
 E  fluxo_espera(2,2,4)
 E  fluxo_espera(2,2,5)
 E  fluxo_espera(2,2,6)
 E  fluxo_espera(2,2,7)
 E  fluxo_espera(2,2,8)
 E  fluxo_espera(2,2,9)
 E  fluxo_espera(2,2,10)
 E  fluxo_espera(3,3,7)
 E  fluxo_espera(3,3,8)
 E  fluxo_espera(3,3,9)
 E  fluxo_espera(3,3,10)
 E  fluxo_espera(3,3,11)
 E  ultimo_tempo(1,1)
 E  ultimo_tempo(1,2)
 E  ultimo_tempo(2,3)
 E  ultimo_tempo(2,1)
 E  ultimo_tempo(3,2)
 E  ultimo_tempo(1,3)
 E  ultimo_tempo(3,1)
 E  ultimo_tempo(2,2)
 E  ultimo_tempo(3,3)
 G  makespan(1)
 G  makespan(2)
 G  makespan(3)
COLUMNS
    MARKER    'MARKER'                 'INTORG'
    f(1,1)    inicio_maquina(m1,t1)  1
    f(1,1)    fluxo_maquina(1,2)  1
    x(1,1,1)  inicio_maquina(m1,t1)  1
    x(1,1,1)  fluxo_maquina(1,2)  1
    x(1,1,1)  inicio_espera(m0,j1,t1)  -1
    x(1,1,1)  fluxo_espera(2,1,2)  1
    e(1,1,1)  inicio_espera(m0,j1,t1)  -1
    e(1,1,1)  fluxo_espera(1,1,2)  1
    x(1,2,1)  inicio_maquina(m1,t1)  1
    x(1,2,1)  fluxo_maquina(1,2)  1
    x(1,2,1)  inicio_espera(m0,j2,t1)  -1
    x(1,2,1)  fluxo_espera(3,2,2)  1
    e(1,2,1)  inicio_espera(m0,j2,t1)  -1
    e(1,2,1)  fluxo_espera(1,2,2)  1
    f(1,2)    fluxo_maquina(1,2)  -1
    f(1,2)    fluxo_maquina(1,3)  1
    x(1,1,2)  fluxo_maquina(1,2)  -1
    x(1,1,2)  fluxo_maquina(1,3)  1
    x(1,1,2)  fluxo_espera(1,1,2)  -1
    x(1,1,2)  fluxo_espera(2,1,3)  1
    e(1,1,2)  fluxo_espera(1,1,2)  -1
    e(1,1,2)  fluxo_espera(1,1,3)  1
    x(1,2,2)  fluxo_maquina(1,2)  -1
    x(1,2,2)  fluxo_maquina(1,3)  1
    x(1,2,2)  fluxo_espera(1,2,2)  -1
    x(1,2,2)  fluxo_espera(3,2,3)  1
    e(1,2,2)  fluxo_espera(1,2,2)  -1
    e(1,2,2)  fluxo_espera(1,2,3)  1
    f(1,3)    fluxo_maquina(1,3)  -1
    f(1,3)    fluxo_maquina(1,4)  1
    x(1,1,3)  fluxo_maquina(1,3)  -1
    x(1,1,3)  fluxo_maquina(1,4)  1
    x(1,1,3)  fluxo_espera(1,1,3)  -1
    x(1,1,3)  fluxo_espera(2,1,4)  1
    e(1,1,3)  fluxo_espera(1,1,3)  -1
    e(1,1,3)  fluxo_espera(1,1,4)  1
    x(1,2,3)  fluxo_maquina(1,3)  -1
    x(1,2,3)  fluxo_maquina(1,4)  1
    x(1,2,3)  fluxo_espera(1,2,3)  -1
    x(1,2,3)  fluxo_espera(3,2,4)  1
    e(1,2,3)  fluxo_espera(1,2,3)  -1
    e(1,2,3)  fluxo_espera(1,2,4)  1
    f(1,4)    fluxo_maquina(1,4)  -1
    f(1,4)    fluxo_maquina(1,5)  1
    x(1,1,4)  fluxo_maquina(1,4)  -1
    x(1,1,4)  fluxo_maquina(1,5)  1
    x(1,1,4)  fluxo_espera(1,1,4)  -1
    x(1,1,4)  fluxo_espera(2,1,5)  1
    e(1,1,4)  fluxo_espera(1,1,4)  -1
    e(1,1,4)  fluxo_espera(1,1,5)  1
    x(1,2,4)  fluxo_maquina(1,4)  -1
    x(1,2,4)  fluxo_maquina(1,5)  1
    x(1,2,4)  fluxo_espera(1,2,4)  -1
    x(1,2,4)  fluxo_espera(3,2,5)  1
    e(1,2,4)  fluxo_espera(1,2,4)  -1
    e(1,2,4)  fluxo_espera(1,2,5)  1
    f(1,5)    fluxo_maquina(1,5)  -1
    f(1,5)    fluxo_maquina(1,6)  1
    x(1,1,5)  fluxo_maquina(1,5)  -1
    x(1,1,5)  fluxo_maquina(1,6)  1
    x(1,1,5)  fluxo_espera(1,1,5)  -1
    x(1,1,5)  fluxo_espera(2,1,6)  1
    e(1,1,5)  fluxo_espera(1,1,5)  -1
    e(1,1,5)  fluxo_espera(1,1,6)  1
    x(1,2,5)  fluxo_maquina(1,5)  -1
    x(1,2,5)  fluxo_maquina(1,6)  1
    x(1,2,5)  fluxo_espera(1,2,5)  -1
    x(1,2,5)  fluxo_espera(3,2,6)  1
    e(1,2,5)  fluxo_espera(1,2,5)  -1
    e(1,2,5)  fluxo_espera(1,2,6)  1
    x(1,3,5)  fluxo_maquina(1,5)  -1
    x(1,3,5)  fluxo_maquina(1,7)  1
    x(1,3,5)  fluxo_espera(1,3,5)  -1
    x(1,3,5)  fluxo_espera(3,3,7)  1
    e(1,3,5)  fluxo_espera(1,3,5)  -1
    e(1,3,5)  fluxo_espera(1,3,6)  1
    f(1,6)    fluxo_maquina(1,6)  -1
    f(1,6)    fluxo_maquina(1,7)  1
    x(1,1,6)  fluxo_maquina(1,6)  -1
    x(1,1,6)  fluxo_maquina(1,7)  1
    x(1,1,6)  fluxo_espera(1,1,6)  -1
    x(1,1,6)  fluxo_espera(2,1,7)  1
    e(1,1,6)  fluxo_espera(1,1,6)  -1
    e(1,1,6)  fluxo_espera(1,1,7)  1
    x(1,2,6)  fluxo_maquina(1,6)  -1
    x(1,2,6)  fluxo_maquina(1,7)  1
    x(1,2,6)  fluxo_espera(1,2,6)  -1
    x(1,2,6)  fluxo_espera(3,2,7)  1
    e(1,2,6)  fluxo_espera(1,2,6)  -1
    e(1,2,6)  fluxo_espera(1,2,7)  1
    x(1,3,6)  fluxo_maquina(1,6)  -1
    x(1,3,6)  fluxo_maquina(1,8)  1
    x(1,3,6)  fluxo_espera(1,3,6)  -1
    x(1,3,6)  fluxo_espera(3,3,8)  1
    e(1,3,6)  fluxo_espera(1,3,6)  -1
    e(1,3,6)  fluxo_espera(1,3,7)  1
    f(1,7)    fluxo_maquina(1,7)  -1
    f(1,7)    fluxo_maquina(1,8)  1
    x(1,1,7)  fluxo_maquina(1,7)  -1
    x(1,1,7)  fluxo_maquina(1,8)  1
    x(1,1,7)  fluxo_espera(1,1,7)  -1
    x(1,1,7)  fluxo_espera(2,1,8)  1
    e(1,1,7)  fluxo_espera(1,1,7)  -1
    e(1,1,7)  ultimo_tempo(1,1)  1
    x(1,2,7)  fluxo_maquina(1,7)  -1
    x(1,2,7)  fluxo_maquina(1,8)  1
    x(1,2,7)  fluxo_espera(1,2,7)  -1
    x(1,2,7)  fluxo_espera(3,2,8)  1
    e(1,2,7)  fluxo_espera(1,2,7)  -1
    e(1,2,7)  ultimo_tempo(1,2)  1
    x(1,3,7)  fluxo_maquina(1,7)  -1
    x(1,3,7)  fluxo_maquina(1,9)  1
    x(1,3,7)  fluxo_espera(1,3,7)  -1
    x(1,3,7)  fluxo_espera(3,3,9)  1
    e(1,3,7)  fluxo_espera(1,3,7)  -1
    e(1,3,7)  fluxo_espera(1,3,8)  1
    f(1,8)    fluxo_maquina(1,8)  -1
    f(1,8)    fluxo_maquina(1,9)  1
    x(1,1,8)  fluxo_maquina(1,8)  -1
    x(1,1,8)  fluxo_maquina(1,9)  1
    x(1,1,8)  ultimo_tempo(1,1)  -1
    x(1,1,8)  ultimo_tempo(2,1)  1
    e(1,1,8)  OBJ       0
    x(1,2,8)  fluxo_maquina(1,8)  -1
    x(1,2,8)  fluxo_maquina(1,9)  1
    x(1,2,8)  ultimo_tempo(1,2)  -1
    x(1,2,8)  ultimo_tempo(3,2)  1
    e(1,2,8)  OBJ       0
    x(1,3,8)  fluxo_maquina(1,8)  -1
    x(1,3,8)  fluxo_maquina(1,10)  1
    x(1,3,8)  fluxo_espera(1,3,8)  -1
    x(1,3,8)  fluxo_espera(3,3,10)  1
    e(1,3,8)  fluxo_espera(1,3,8)  -1
    e(1,3,8)  fluxo_espera(1,3,9)  1
    f(1,9)    fluxo_maquina(1,9)  -1
    f(1,9)    fluxo_maquina(1,10)  1
    x(1,3,9)  fluxo_maquina(1,9)  -1
    x(1,3,9)  fluxo_maquina(1,11)  1
    x(1,3,9)  fluxo_espera(1,3,9)  -1
    x(1,3,9)  fluxo_espera(3,3,11)  1
    e(1,3,9)  fluxo_espera(1,3,9)  -1
    e(1,3,9)  ultimo_tempo(1,3)  1
    f(1,10)   fluxo_maquina(1,10)  -1
    f(1,10)   fluxo_maquina(1,11)  1
    x(1,3,10)  fluxo_maquina(1,10)  -1
    x(1,3,10)  fluxo_maquina(1,12)  1
    x(1,3,10)  ultimo_tempo(1,3)  -1
    x(1,3,10)  ultimo_tempo(3,3)  1
    e(1,3,10)  OBJ       0
    f(1,11)   fluxo_maquina(1,11)  -1
    f(1,11)   fluxo_maquina(1,12)  1
    f(1,12)   fluxo_maquina(1,12)  -1
    f(2,1)    inicio_maquina(m2,t1)  1
    f(2,1)    fluxo_maquina(2,2)  1
    x(2,3,1)  inicio_maquina(m2,t1)  1
    x(2,3,1)  fluxo_maquina(2,5)  1
    x(2,3,1)  inicio_espera(m1,j3,t1)  -1
    x(2,3,1)  fluxo_espera(1,3,5)  1
    e(2,3,1)  inicio_espera(m1,j3,t1)  -1
    e(2,3,1)  fluxo_espera(2,3,2)  1
    f(2,2)    fluxo_maquina(2,2)  -1
    f(2,2)    fluxo_maquina(2,3)  1
    x(2,1,2)  fluxo_maquina(2,2)  -1
    x(2,1,2)  fluxo_maquina(2,4)  1
    x(2,1,2)  fluxo_espera(2,1,2)  -1
    x(2,1,2)  fluxo_espera(3,1,4)  1
    e(2,1,2)  fluxo_espera(2,1,2)  -1
    e(2,1,2)  fluxo_espera(2,1,3)  1
    x(2,3,2)  fluxo_maquina(2,2)  -1
    x(2,3,2)  fluxo_maquina(2,6)  1
    x(2,3,2)  fluxo_espera(2,3,2)  -1
    x(2,3,2)  fluxo_espera(1,3,6)  1
    e(2,3,2)  fluxo_espera(2,3,2)  -1
    e(2,3,2)  fluxo_espera(2,3,3)  1
    f(2,3)    fluxo_maquina(2,3)  -1
    f(2,3)    fluxo_maquina(2,4)  1
    x(2,1,3)  fluxo_maquina(2,3)  -1
    x(2,1,3)  fluxo_maquina(2,5)  1
    x(2,1,3)  fluxo_espera(2,1,3)  -1
    x(2,1,3)  fluxo_espera(3,1,5)  1
    e(2,1,3)  fluxo_espera(2,1,3)  -1
    e(2,1,3)  fluxo_espera(2,1,4)  1
    x(2,3,3)  fluxo_maquina(2,3)  -1
    x(2,3,3)  fluxo_maquina(2,7)  1
    x(2,3,3)  fluxo_espera(2,3,3)  -1
    x(2,3,3)  fluxo_espera(1,3,7)  1
    e(2,3,3)  fluxo_espera(2,3,3)  -1
    e(2,3,3)  fluxo_espera(2,3,4)  1
    f(2,4)    fluxo_maquina(2,4)  -1
    f(2,4)    fluxo_maquina(2,5)  1
    x(2,1,4)  fluxo_maquina(2,4)  -1
    x(2,1,4)  fluxo_maquina(2,6)  1
    x(2,1,4)  fluxo_espera(2,1,4)  -1
    x(2,1,4)  fluxo_espera(3,1,6)  1
    e(2,1,4)  fluxo_espera(2,1,4)  -1
    e(2,1,4)  fluxo_espera(2,1,5)  1
    x(2,2,4)  fluxo_maquina(2,4)  -1
    x(2,2,4)  fluxo_maquina(2,6)  1
    x(2,2,4)  fluxo_espera(2,2,4)  -1
    x(2,2,4)  makespan(2)  -5
    e(2,2,4)  fluxo_espera(2,2,4)  -1
    e(2,2,4)  fluxo_espera(2,2,5)  1
    x(2,3,4)  fluxo_maquina(2,4)  -1
    x(2,3,4)  fluxo_maquina(2,8)  1
    x(2,3,4)  fluxo_espera(2,3,4)  -1
    x(2,3,4)  fluxo_espera(1,3,8)  1
    e(2,3,4)  fluxo_espera(2,3,4)  -1
    e(2,3,4)  fluxo_espera(2,3,5)  1
    f(2,5)    fluxo_maquina(2,5)  -1
    f(2,5)    fluxo_maquina(2,6)  1
    x(2,1,5)  fluxo_maquina(2,5)  -1
    x(2,1,5)  fluxo_maquina(2,7)  1
    x(2,1,5)  fluxo_espera(2,1,5)  -1
    x(2,1,5)  fluxo_espera(3,1,7)  1
    e(2,1,5)  fluxo_espera(2,1,5)  -1
    e(2,1,5)  fluxo_espera(2,1,6)  1
    x(2,2,5)  fluxo_maquina(2,5)  -1
    x(2,2,5)  fluxo_maquina(2,7)  1
    x(2,2,5)  fluxo_espera(2,2,5)  -1
    x(2,2,5)  makespan(2)  -6
    e(2,2,5)  fluxo_espera(2,2,5)  -1
    e(2,2,5)  fluxo_espera(2,2,6)  1
    x(2,3,5)  fluxo_maquina(2,5)  -1
    x(2,3,5)  fluxo_maquina(2,9)  1
    x(2,3,5)  fluxo_espera(2,3,5)  -1
    x(2,3,5)  fluxo_espera(1,3,9)  1
    e(2,3,5)  fluxo_espera(2,3,5)  -1
    e(2,3,5)  ultimo_tempo(2,3)  1
    f(2,6)    fluxo_maquina(2,6)  -1
    f(2,6)    fluxo_maquina(2,7)  1
    x(2,1,6)  fluxo_maquina(2,6)  -1
    x(2,1,6)  fluxo_maquina(2,8)  1
    x(2,1,6)  fluxo_espera(2,1,6)  -1
    x(2,1,6)  fluxo_espera(3,1,8)  1
    e(2,1,6)  fluxo_espera(2,1,6)  -1
    e(2,1,6)  fluxo_espera(2,1,7)  1
    x(2,2,6)  fluxo_maquina(2,6)  -1
    x(2,2,6)  fluxo_maquina(2,8)  1
    x(2,2,6)  fluxo_espera(2,2,6)  -1
    x(2,2,6)  makespan(2)  -7
    e(2,2,6)  fluxo_espera(2,2,6)  -1
    e(2,2,6)  fluxo_espera(2,2,7)  1
    x(2,3,6)  fluxo_maquina(2,6)  -1
    x(2,3,6)  fluxo_maquina(2,10)  1
    x(2,3,6)  ultimo_tempo(2,3)  -1
    x(2,3,6)  ultimo_tempo(1,3)  1
    e(2,3,6)  OBJ       0
    f(2,7)    fluxo_maquina(2,7)  -1
    f(2,7)    fluxo_maquina(2,8)  1
    x(2,1,7)  fluxo_maquina(2,7)  -1
    x(2,1,7)  fluxo_maquina(2,9)  1
    x(2,1,7)  fluxo_espera(2,1,7)  -1
    x(2,1,7)  fluxo_espera(3,1,9)  1
    e(2,1,7)  fluxo_espera(2,1,7)  -1
    e(2,1,7)  fluxo_espera(2,1,8)  1
    x(2,2,7)  fluxo_maquina(2,7)  -1
    x(2,2,7)  fluxo_maquina(2,9)  1
    x(2,2,7)  fluxo_espera(2,2,7)  -1
    x(2,2,7)  makespan(2)  -8
    e(2,2,7)  fluxo_espera(2,2,7)  -1
    e(2,2,7)  fluxo_espera(2,2,8)  1
    f(2,8)    fluxo_maquina(2,8)  -1
    f(2,8)    fluxo_maquina(2,9)  1
    x(2,1,8)  fluxo_maquina(2,8)  -1
    x(2,1,8)  fluxo_maquina(2,10)  1
    x(2,1,8)  fluxo_espera(2,1,8)  -1
    x(2,1,8)  fluxo_espera(3,1,10)  1
    e(2,1,8)  fluxo_espera(2,1,8)  -1
    e(2,1,8)  ultimo_tempo(2,1)  1
    x(2,2,8)  fluxo_maquina(2,8)  -1
    x(2,2,8)  fluxo_maquina(2,10)  1
    x(2,2,8)  fluxo_espera(2,2,8)  -1
    x(2,2,8)  makespan(2)  -9
    e(2,2,8)  fluxo_espera(2,2,8)  -1
    e(2,2,8)  fluxo_espera(2,2,9)  1
    f(2,9)    fluxo_maquina(2,9)  -1
    f(2,9)    fluxo_maquina(2,10)  1
    x(2,1,9)  fluxo_maquina(2,9)  -1
    x(2,1,9)  fluxo_maquina(2,11)  1
    x(2,1,9)  ultimo_tempo(2,1)  -1
    x(2,1,9)  ultimo_tempo(3,1)  1
    e(2,1,9)  OBJ       0
    x(2,2,9)  fluxo_maquina(2,9)  -1
    x(2,2,9)  fluxo_maquina(2,11)  1
    x(2,2,9)  fluxo_espera(2,2,9)  -1
    x(2,2,9)  makespan(2)  -10
    e(2,2,9)  fluxo_espera(2,2,9)  -1
    e(2,2,9)  fluxo_espera(2,2,10)  1
    f(2,10)   fluxo_maquina(2,10)  -1
    f(2,10)   fluxo_maquina(2,11)  1
    x(2,2,10)  fluxo_maquina(2,10)  -1
    x(2,2,10)  fluxo_maquina(2,12)  1
    x(2,2,10)  fluxo_espera(2,2,10)  -1
    x(2,2,10)  makespan(2)  -11
    e(2,2,10)  fluxo_espera(2,2,10)  -1
    e(2,2,10)  ultimo_tempo(2,2)  1
    f(2,11)   fluxo_maquina(2,11)  -1
    f(2,11)   fluxo_maquina(2,12)  1
    x(2,2,11)  fluxo_maquina(2,11)  -1
    x(2,2,11)  ultimo_tempo(2,2)  -1
    x(2,2,11)  makespan(2)  -12
    e(2,2,11)  OBJ       0
    f(2,12)   fluxo_maquina(2,12)  -1
    f(3,1)    inicio_maquina(m3,t1)  1
    f(3,1)    fluxo_maquina(3,2)  1
    f(3,2)    fluxo_maquina(3,2)  -1
    f(3,2)    fluxo_maquina(3,3)  1
    x(3,2,2)  fluxo_maquina(3,2)  -1
    x(3,2,2)  fluxo_maquina(3,4)  1
    x(3,2,2)  fluxo_espera(3,2,2)  -1
    x(3,2,2)  fluxo_espera(2,2,4)  1
    e(3,2,2)  fluxo_espera(3,2,2)  -1
    e(3,2,2)  fluxo_espera(3,2,3)  1
    f(3,3)    fluxo_maquina(3,3)  -1
    f(3,3)    fluxo_maquina(3,4)  1
    x(3,2,3)  fluxo_maquina(3,3)  -1
    x(3,2,3)  fluxo_maquina(3,5)  1
    x(3,2,3)  fluxo_espera(3,2,3)  -1
    x(3,2,3)  fluxo_espera(2,2,5)  1
    e(3,2,3)  fluxo_espera(3,2,3)  -1
    e(3,2,3)  fluxo_espera(3,2,4)  1
    f(3,4)    fluxo_maquina(3,4)  -1
    f(3,4)    fluxo_maquina(3,5)  1
    x(3,1,4)  fluxo_maquina(3,4)  -1
    x(3,1,4)  fluxo_maquina(3,6)  1
    x(3,1,4)  fluxo_espera(3,1,4)  -1
    x(3,1,4)  makespan(1)  -5
    e(3,1,4)  fluxo_espera(3,1,4)  -1
    e(3,1,4)  fluxo_espera(3,1,5)  1
    x(3,2,4)  fluxo_maquina(3,4)  -1
    x(3,2,4)  fluxo_maquina(3,6)  1
    x(3,2,4)  fluxo_espera(3,2,4)  -1
    x(3,2,4)  fluxo_espera(2,2,6)  1
    e(3,2,4)  fluxo_espera(3,2,4)  -1
    e(3,2,4)  fluxo_espera(3,2,5)  1
    f(3,5)    fluxo_maquina(3,5)  -1
    f(3,5)    fluxo_maquina(3,6)  1
    x(3,1,5)  fluxo_maquina(3,5)  -1
    x(3,1,5)  fluxo_maquina(3,7)  1
    x(3,1,5)  fluxo_espera(3,1,5)  -1
    x(3,1,5)  makespan(1)  -6
    e(3,1,5)  fluxo_espera(3,1,5)  -1
    e(3,1,5)  fluxo_espera(3,1,6)  1
    x(3,2,5)  fluxo_maquina(3,5)  -1
    x(3,2,5)  fluxo_maquina(3,7)  1
    x(3,2,5)  fluxo_espera(3,2,5)  -1
    x(3,2,5)  fluxo_espera(2,2,7)  1
    e(3,2,5)  fluxo_espera(3,2,5)  -1
    e(3,2,5)  fluxo_espera(3,2,6)  1
    f(3,6)    fluxo_maquina(3,6)  -1
    f(3,6)    fluxo_maquina(3,7)  1
    x(3,1,6)  fluxo_maquina(3,6)  -1
    x(3,1,6)  fluxo_maquina(3,8)  1
    x(3,1,6)  fluxo_espera(3,1,6)  -1
    x(3,1,6)  makespan(1)  -7
    e(3,1,6)  fluxo_espera(3,1,6)  -1
    e(3,1,6)  fluxo_espera(3,1,7)  1
    x(3,2,6)  fluxo_maquina(3,6)  -1
    x(3,2,6)  fluxo_maquina(3,8)  1
    x(3,2,6)  fluxo_espera(3,2,6)  -1
    x(3,2,6)  fluxo_espera(2,2,8)  1
    e(3,2,6)  fluxo_espera(3,2,6)  -1
    e(3,2,6)  fluxo_espera(3,2,7)  1
    f(3,7)    fluxo_maquina(3,7)  -1
    f(3,7)    fluxo_maquina(3,8)  1
    x(3,1,7)  fluxo_maquina(3,7)  -1
    x(3,1,7)  fluxo_maquina(3,9)  1
    x(3,1,7)  fluxo_espera(3,1,7)  -1
    x(3,1,7)  makespan(1)  -8
    e(3,1,7)  fluxo_espera(3,1,7)  -1
    e(3,1,7)  fluxo_espera(3,1,8)  1
    x(3,2,7)  fluxo_maquina(3,7)  -1
    x(3,2,7)  fluxo_maquina(3,9)  1
    x(3,2,7)  fluxo_espera(3,2,7)  -1
    x(3,2,7)  fluxo_espera(2,2,9)  1
    e(3,2,7)  fluxo_espera(3,2,7)  -1
    e(3,2,7)  fluxo_espera(3,2,8)  1
    x(3,3,7)  fluxo_maquina(3,7)  -1
    x(3,3,7)  fluxo_maquina(3,8)  1
    x(3,3,7)  fluxo_espera(3,3,7)  -1
    x(3,3,7)  makespan(3)  -7
    e(3,3,7)  fluxo_espera(3,3,7)  -1
    e(3,3,7)  fluxo_espera(3,3,8)  1
    f(3,8)    fluxo_maquina(3,8)  -1
    f(3,8)    fluxo_maquina(3,9)  1
    x(3,1,8)  fluxo_maquina(3,8)  -1
    x(3,1,8)  fluxo_maquina(3,10)  1
    x(3,1,8)  fluxo_espera(3,1,8)  -1
    x(3,1,8)  makespan(1)  -9
    e(3,1,8)  fluxo_espera(3,1,8)  -1
    e(3,1,8)  fluxo_espera(3,1,9)  1
    x(3,2,8)  fluxo_maquina(3,8)  -1
    x(3,2,8)  fluxo_maquina(3,10)  1
    x(3,2,8)  fluxo_espera(3,2,8)  -1
    x(3,2,8)  fluxo_espera(2,2,10)  1
    e(3,2,8)  fluxo_espera(3,2,8)  -1
    e(3,2,8)  ultimo_tempo(3,2)  1
    x(3,3,8)  fluxo_maquina(3,8)  -1
    x(3,3,8)  fluxo_maquina(3,9)  1
    x(3,3,8)  fluxo_espera(3,3,8)  -1
    x(3,3,8)  makespan(3)  -8
    e(3,3,8)  fluxo_espera(3,3,8)  -1
    e(3,3,8)  fluxo_espera(3,3,9)  1
    f(3,9)    fluxo_maquina(3,9)  -1
    f(3,9)    fluxo_maquina(3,10)  1
    x(3,1,9)  fluxo_maquina(3,9)  -1
    x(3,1,9)  fluxo_maquina(3,11)  1
    x(3,1,9)  fluxo_espera(3,1,9)  -1
    x(3,1,9)  makespan(1)  -10
    e(3,1,9)  fluxo_espera(3,1,9)  -1
    e(3,1,9)  fluxo_espera(3,1,10)  1
    x(3,2,9)  fluxo_maquina(3,9)  -1
    x(3,2,9)  fluxo_maquina(3,11)  1
    x(3,2,9)  ultimo_tempo(3,2)  -1
    x(3,2,9)  ultimo_tempo(2,2)  1
    e(3,2,9)  OBJ       0
    x(3,3,9)  fluxo_maquina(3,9)  -1
    x(3,3,9)  fluxo_maquina(3,10)  1
    x(3,3,9)  fluxo_espera(3,3,9)  -1
    x(3,3,9)  makespan(3)  -9
    e(3,3,9)  fluxo_espera(3,3,9)  -1
    e(3,3,9)  fluxo_espera(3,3,10)  1
    f(3,10)   fluxo_maquina(3,10)  -1
    f(3,10)   fluxo_maquina(3,11)  1
    x(3,1,10)  fluxo_maquina(3,10)  -1
    x(3,1,10)  fluxo_maquina(3,12)  1
    x(3,1,10)  fluxo_espera(3,1,10)  -1
    x(3,1,10)  makespan(1)  -11
    e(3,1,10)  fluxo_espera(3,1,10)  -1
    e(3,1,10)  ultimo_tempo(3,1)  1
    x(3,3,10)  fluxo_maquina(3,10)  -1
    x(3,3,10)  fluxo_maquina(3,11)  1
    x(3,3,10)  fluxo_espera(3,3,10)  -1
    x(3,3,10)  makespan(3)  -10
    e(3,3,10)  fluxo_espera(3,3,10)  -1
    e(3,3,10)  fluxo_espera(3,3,11)  1
    f(3,11)   fluxo_maquina(3,11)  -1
    f(3,11)   fluxo_maquina(3,12)  1
    x(3,1,11)  fluxo_maquina(3,11)  -1
    x(3,1,11)  ultimo_tempo(3,1)  -1
    x(3,1,11)  makespan(1)  -12
    e(3,1,11)  OBJ       0
    x(3,3,11)  fluxo_maquina(3,11)  -1
    x(3,3,11)  fluxo_maquina(3,12)  1
    x(3,3,11)  fluxo_espera(3,3,11)  -1
    x(3,3,11)  makespan(3)  -11
    e(3,3,11)  fluxo_espera(3,3,11)  -1
    e(3,3,11)  ultimo_tempo(3,3)  1
    f(3,12)   fluxo_maquina(3,12)  -1
    x(3,3,12)  fluxo_maquina(3,12)  -1
    x(3,3,12)  ultimo_tempo(3,3)  -1
    x(3,3,12)  makespan(3)  -12
    e(3,3,12)  OBJ       0
    C         OBJ       1
    C         makespan(1)  1
    C         makespan(2)  1
    C         makespan(3)  1
    MARKER    'MARKER'                 'INTEND'
RHS
    RHS1      inicio_maquina(m1,t1)  1
    RHS1      inicio_maquina(m2,t1)  1
    RHS1      inicio_maquina(m3,t1)  1
    RHS1      inicio_espera(m0,j1,t1)  -1
    RHS1      inicio_espera(m0,j2,t1)  -1
    RHS1      inicio_espera(m1,j3,t1)  -1
BOUNDS
 BV BND1      f(1,1)  
 BV BND1      x(1,1,1)
 BV BND1      e(1,1,1)
 BV BND1      x(1,2,1)
 BV BND1      e(1,2,1)
 BV BND1      f(1,2)  
 BV BND1      x(1,1,2)
 BV BND1      e(1,1,2)
 BV BND1      x(1,2,2)
 BV BND1      e(1,2,2)
 BV BND1      f(1,3)  
 BV BND1      x(1,1,3)
 BV BND1      e(1,1,3)
 BV BND1      x(1,2,3)
 BV BND1      e(1,2,3)
 BV BND1      f(1,4)  
 BV BND1      x(1,1,4)
 BV BND1      e(1,1,4)
 BV BND1      x(1,2,4)
 BV BND1      e(1,2,4)
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
 BV BND1      x(1,3,9)
 BV BND1      e(1,3,9)
 BV BND1      f(1,10) 
 BV BND1      x(1,3,10)
 BV BND1      e(1,3,10)
 BV BND1      f(1,11) 
 BV BND1      f(1,12) 
 BV BND1      f(2,1)  
 BV BND1      x(2,3,1)
 BV BND1      e(2,3,1)
 BV BND1      f(2,2)  
 BV BND1      x(2,1,2)
 BV BND1      e(2,1,2)
 BV BND1      x(2,3,2)
 BV BND1      e(2,3,2)
 BV BND1      f(2,3)  
 BV BND1      x(2,1,3)
 BV BND1      e(2,1,3)
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
 BV BND1      x(2,2,7)
 BV BND1      e(2,2,7)
 BV BND1      f(2,8)  
 BV BND1      x(2,1,8)
 BV BND1      e(2,1,8)
 BV BND1      x(2,2,8)
 BV BND1      e(2,2,8)
 BV BND1      f(2,9)  
 BV BND1      x(2,1,9)
 BV BND1      e(2,1,9)
 BV BND1      x(2,2,9)
 BV BND1      e(2,2,9)
 BV BND1      f(2,10) 
 BV BND1      x(2,2,10)
 BV BND1      e(2,2,10)
 BV BND1      f(2,11) 
 BV BND1      x(2,2,11)
 BV BND1      e(2,2,11)
 BV BND1      f(2,12) 
 BV BND1      f(3,1)  
 BV BND1      f(3,2)  
 BV BND1      x(3,2,2)
 BV BND1      e(3,2,2)
 BV BND1      f(3,3)  
 BV BND1      x(3,2,3)
 BV BND1      e(3,2,3)
 BV BND1      f(3,4)  
 BV BND1      x(3,1,4)
 BV BND1      e(3,1,4)
 BV BND1      x(3,2,4)
 BV BND1      e(3,2,4)
 BV BND1      f(3,5)  
 BV BND1      x(3,1,5)
 BV BND1      e(3,1,5)
 BV BND1      x(3,2,5)
 BV BND1      e(3,2,5)
 BV BND1      f(3,6)  
 BV BND1      x(3,1,6)
 BV BND1      e(3,1,6)
 BV BND1      x(3,2,6)
 BV BND1      e(3,2,6)
 BV BND1      f(3,7)  
 BV BND1      x(3,1,7)
 BV BND1      e(3,1,7)
 BV BND1      x(3,2,7)
 BV BND1      e(3,2,7)
 BV BND1      x(3,3,7)
 BV BND1      e(3,3,7)
 BV BND1      f(3,8)  
 BV BND1      x(3,1,8)
 BV BND1      e(3,1,8)
 BV BND1      x(3,2,8)
 BV BND1      e(3,2,8)
 BV BND1      x(3,3,8)
 BV BND1      e(3,3,8)
 BV BND1      f(3,9)  
 BV BND1      x(3,1,9)
 BV BND1      e(3,1,9)
 BV BND1      x(3,2,9)
 BV BND1      e(3,2,9)
 BV BND1      x(3,3,9)
 BV BND1      e(3,3,9)
 BV BND1      f(3,10) 
 BV BND1      x(3,1,10)
 BV BND1      e(3,1,10)
 BV BND1      x(3,3,10)
 BV BND1      e(3,3,10)
 BV BND1      f(3,11) 
 BV BND1      x(3,1,11)
 BV BND1      e(3,1,11)
 BV BND1      x(3,3,11)
 BV BND1      e(3,3,11)
 BV BND1      f(3,12) 
 BV BND1      x(3,3,12)
 BV BND1      e(3,3,12)
 LI BND1      C         0
ENDATA
