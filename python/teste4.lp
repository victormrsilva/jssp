\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  x(j0,m1) - x(j0,m2) >= 4
ord(0,2):  x(j0,m0) - x(j0,m1) >= 4
ord(1,1):  - x(j1,m0) + x(j1,m1) >= 10
ord(1,2):  - x(j1,m1) + x(j1,m2) >= 3
ord(2,1):  x(j2,m0) - x(j2,m1) >= 4
ord(2,2):  - x(j2,m0) + x(j2,m2) >= 1
phi(0,1,0):  x(j0,m0) + 38 y(j0,k1,m0) - x(j1,m0) >= 10
psy(0,1,0):  - x(j0,m0) -38 y(j0,k1,m0) + x(j1,m0) >= -37
phi(0,1,1):  x(j0,m1) + 38 y(j0,k1,m1) - x(j1,m1) >= 3
psy(0,1,1):  - x(j0,m1) -38 y(j0,k1,m1) + x(j1,m1) >= -34
phi(0,1,2):  x(j0,m2) + 38 y(j0,k1,m2) - x(j1,m2) >= 5
psy(0,1,2):  - x(j0,m2) -38 y(j0,k1,m2) + x(j1,m2) >= -34
phi(0,2,0):  x(j0,m0) + 38 y(j0,k2,m0) - x(j2,m0) >= 1
psy(0,2,0):  - x(j0,m0) -38 y(j0,k2,m0) + x(j2,m0) >= -37
phi(0,2,1):  x(j0,m1) + 38 y(j0,k2,m1) - x(j2,m1) >= 4
psy(0,2,1):  - x(j0,m1) -38 y(j0,k2,m1) + x(j2,m1) >= -34
phi(0,2,2):  x(j0,m2) + 38 y(j0,k2,m2) - x(j2,m2) >= 6
psy(0,2,2):  - x(j0,m2) -38 y(j0,k2,m2) + x(j2,m2) >= -34
phi(1,0,0):  - x(j0,m0) + 38 y(j1,k0,m0) + x(j1,m0) >= 1
psy(1,0,0):  x(j0,m0) -38 y(j1,k0,m0) - x(j1,m0) >= -28
phi(1,0,1):  - x(j0,m1) + 38 y(j1,k0,m1) + x(j1,m1) >= 4
psy(1,0,1):  x(j0,m1) -38 y(j1,k0,m1) - x(j1,m1) >= -35
phi(1,0,2):  - x(j0,m2) + 38 y(j1,k0,m2) + x(j1,m2) >= 4
psy(1,0,2):  x(j0,m2) -38 y(j1,k0,m2) - x(j1,m2) >= -33
phi(1,2,0):  x(j1,m0) + 38 y(j1,k2,m0) - x(j2,m0) >= 1
psy(1,2,0):  - x(j1,m0) -38 y(j1,k2,m0) + x(j2,m0) >= -28
phi(1,2,1):  x(j1,m1) + 38 y(j1,k2,m1) - x(j2,m1) >= 4
psy(1,2,1):  - x(j1,m1) -38 y(j1,k2,m1) + x(j2,m1) >= -35
phi(1,2,2):  x(j1,m2) + 38 y(j1,k2,m2) - x(j2,m2) >= 6
psy(1,2,2):  - x(j1,m2) -38 y(j1,k2,m2) + x(j2,m2) >= -33
phi(2,0,0):  - x(j0,m0) + 38 y(j2,k0,m0) + x(j2,m0) >= 1
psy(2,0,0):  x(j0,m0) -38 y(j2,k0,m0) - x(j2,m0) >= -37
phi(2,0,1):  - x(j0,m1) + 38 y(j2,k0,m1) + x(j2,m1) >= 4
psy(2,0,1):  x(j0,m1) -38 y(j2,k0,m1) - x(j2,m1) >= -34
phi(2,0,2):  - x(j0,m2) + 38 y(j2,k0,m2) + x(j2,m2) >= 4
psy(2,0,2):  x(j0,m2) -38 y(j2,k0,m2) - x(j2,m2) >= -32
phi(2,1,0):  - x(j1,m0) + 38 y(j2,k1,m0) + x(j2,m0) >= 10
psy(2,1,0):  x(j1,m0) -38 y(j2,k1,m0) - x(j2,m0) >= -37
phi(2,1,1):  - x(j1,m1) + 38 y(j2,k1,m1) + x(j2,m1) >= 3
psy(2,1,1):  x(j1,m1) -38 y(j2,k1,m1) - x(j2,m1) >= -34
phi(2,1,2):  - x(j1,m2) + 38 y(j2,k1,m2) + x(j2,m2) >= 5
psy(2,1,2):  x(j1,m2) -38 y(j2,k1,m2) - x(j2,m2) >= -32
makespan(0):  C - x(j0,m0) >= 1
makespan(1):  C - x(j1,m2) >= 5
makespan(2):  C - x(j2,m2) >= 6
Bounds
 0 <= C <= 38
 8 <= x(j0,m0) <= 37
 0 <= y(j0,k1,m0) <= 1
 0 <= y(j1,k0,m0) <= 1
 0 <= y(j0,k2,m0) <= 1
 0 <= y(j2,k0,m0) <= 1
 0 <= x(j1,m0) <= 20
 0 <= y(j1,k2,m0) <= 1
 0 <= y(j2,k1,m0) <= 1
 4 <= x(j2,m0) <= 31
 4 <= x(j0,m1) <= 33
 0 <= y(j0,k1,m1) <= 1
 0 <= y(j1,k0,m1) <= 1
 0 <= y(j0,k2,m1) <= 1
 0 <= y(j2,k0,m1) <= 1
 10 <= x(j1,m1) <= 30
 0 <= y(j1,k2,m1) <= 1
 0 <= y(j2,k1,m1) <= 1
 0 <= x(j2,m1) <= 27
 0 <= x(j0,m2) <= 29
 0 <= y(j0,k1,m2) <= 1
 0 <= y(j1,k0,m2) <= 1
 0 <= y(j0,k2,m2) <= 1
 0 <= y(j2,k0,m2) <= 1
 13 <= x(j1,m2) <= 33
 0 <= y(j1,k2,m2) <= 1
 0 <= y(j2,k1,m2) <= 1
 5 <= x(j2,m2) <= 32
Integers
C x(j0,m0) y(j0,k1,m0) y(j1,k0,m0) y(j0,k2,m0) y(j2,k0,m0) x(j1,m0) y(j1,k2,m0) y(j2,k1,m0) x(j2,m0) 
x(j0,m1) y(j0,k1,m1) y(j1,k0,m1) y(j0,k2,m1) y(j2,k0,m1) x(j1,m1) y(j1,k2,m1) y(j2,k1,m1) x(j2,m1) x(j0,m2) 
y(j0,k1,m2) y(j1,k0,m2) y(j0,k2,m2) y(j2,k0,m2) x(j1,m2) y(j1,k2,m2) y(j2,k1,m2) x(j2,m2) 
End
