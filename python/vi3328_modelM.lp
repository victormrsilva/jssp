\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(j0,m1) + x(j0,m2) >= 7
ord(0,2):  x(j0,m0) - x(j0,m2) >= 7
ord(1,1):  x(j1,m0) - x(j1,m1) >= 4
ord(1,2):  - x(j1,m0) + x(j1,m2) >= 3
ord(2,1):  - x(j2,m1) + x(j2,m2) >= 7
ord(2,2):  x(j2,m0) - x(j2,m2) >= 6
phi(0,1,0):  x(j0,m0) + 44 y(j0,k1,m0) - x(j1,m0) >= 3
psy(0,1,0):  - x(j0,m0) -44 y(j0,k1,m0) + x(j1,m0) >= -43
phi(0,1,1):  x(j0,m1) + 44 y(j0,k1,m1) - x(j1,m1) >= 4
psy(0,1,1):  - x(j0,m1) -44 y(j0,k1,m1) + x(j1,m1) >= -37
phi(0,1,2):  x(j0,m2) + 44 y(j0,k1,m2) - x(j1,m2) >= 7
psy(0,1,2):  - x(j0,m2) -44 y(j0,k1,m2) + x(j1,m2) >= -37
phi(0,2,0):  x(j0,m0) + 44 y(j0,k2,m0) - x(j2,m0) >= 2
psy(0,2,0):  - x(j0,m0) -44 y(j0,k2,m0) + x(j2,m0) >= -43
phi(0,2,1):  x(j0,m1) + 44 y(j0,k2,m1) - x(j2,m1) >= 7
psy(0,2,1):  - x(j0,m1) -44 y(j0,k2,m1) + x(j2,m1) >= -37
phi(0,2,2):  x(j0,m2) + 44 y(j0,k2,m2) - x(j2,m2) >= 6
psy(0,2,2):  - x(j0,m2) -44 y(j0,k2,m2) + x(j2,m2) >= -37
phi(1,0,0):  - x(j0,m0) + 44 y(j1,k0,m0) + x(j1,m0) >= 1
psy(1,0,0):  x(j0,m0) -44 y(j1,k0,m0) - x(j1,m0) >= -41
phi(1,0,1):  - x(j0,m1) + 44 y(j1,k0,m1) + x(j1,m1) >= 7
psy(1,0,1):  x(j0,m1) -44 y(j1,k0,m1) - x(j1,m1) >= -40
phi(1,0,2):  - x(j0,m2) + 44 y(j1,k0,m2) + x(j1,m2) >= 7
psy(1,0,2):  x(j0,m2) -44 y(j1,k0,m2) - x(j1,m2) >= -37
phi(1,2,0):  x(j1,m0) + 44 y(j1,k2,m0) - x(j2,m0) >= 2
psy(1,2,0):  - x(j1,m0) -44 y(j1,k2,m0) + x(j2,m0) >= -41
phi(1,2,1):  x(j1,m1) + 44 y(j1,k2,m1) - x(j2,m1) >= 7
psy(1,2,1):  - x(j1,m1) -44 y(j1,k2,m1) + x(j2,m1) >= -40
phi(1,2,2):  x(j1,m2) + 44 y(j1,k2,m2) - x(j2,m2) >= 6
psy(1,2,2):  - x(j1,m2) -44 y(j1,k2,m2) + x(j2,m2) >= -37
phi(2,0,0):  - x(j0,m0) + 44 y(j2,k0,m0) + x(j2,m0) >= 1
psy(2,0,0):  x(j0,m0) -44 y(j2,k0,m0) - x(j2,m0) >= -42
phi(2,0,1):  - x(j0,m1) + 44 y(j2,k0,m1) + x(j2,m1) >= 7
psy(2,0,1):  x(j0,m1) -44 y(j2,k0,m1) - x(j2,m1) >= -37
phi(2,0,2):  - x(j0,m2) + 44 y(j2,k0,m2) + x(j2,m2) >= 7
psy(2,0,2):  x(j0,m2) -44 y(j2,k0,m2) - x(j2,m2) >= -38
phi(2,1,0):  - x(j1,m0) + 44 y(j2,k1,m0) + x(j2,m0) >= 3
psy(2,1,0):  x(j1,m0) -44 y(j2,k1,m0) - x(j2,m0) >= -42
phi(2,1,1):  - x(j1,m1) + 44 y(j2,k1,m1) + x(j2,m1) >= 4
psy(2,1,1):  x(j1,m1) -44 y(j2,k1,m1) - x(j2,m1) >= -37
phi(2,1,2):  - x(j1,m2) + 44 y(j2,k1,m2) + x(j2,m2) >= 7
psy(2,1,2):  x(j1,m2) -44 y(j2,k1,m2) - x(j2,m2) >= -38
makespan(0):  C - x(j0,m0) >= 1
makespan(1):  C - x(j1,m2) >= 7
makespan(2):  C - x(j2,m0) >= 2
Bounds
 0 <= y(j0,k1,m0) <= 1
 0 <= y(j1,k0,m0) <= 1
 0 <= y(j0,k2,m0) <= 1
 0 <= y(j2,k0,m0) <= 1
 0 <= y(j1,k2,m0) <= 1
 0 <= y(j2,k1,m0) <= 1
 0 <= y(j0,k1,m1) <= 1
 0 <= y(j1,k0,m1) <= 1
 0 <= y(j0,k2,m1) <= 1
 0 <= y(j2,k0,m1) <= 1
 0 <= y(j1,k2,m1) <= 1
 0 <= y(j2,k1,m1) <= 1
 0 <= y(j0,k1,m2) <= 1
 0 <= y(j1,k0,m2) <= 1
 0 <= y(j0,k2,m2) <= 1
 0 <= y(j2,k0,m2) <= 1
 0 <= y(j1,k2,m2) <= 1
 0 <= y(j2,k1,m2) <= 1
Integers
C x(j0,m0) y(j0,k1,m0) y(j1,k0,m0) y(j0,k2,m0) y(j2,k0,m0) x(j1,m0) y(j1,k2,m0) y(j2,k1,m0) x(j2,m0) 
x(j0,m1) y(j0,k1,m1) y(j1,k0,m1) y(j0,k2,m1) y(j2,k0,m1) x(j1,m1) y(j1,k2,m1) y(j2,k1,m1) x(j2,m1) x(j0,m2) 
y(j0,k1,m2) y(j1,k0,m2) y(j0,k2,m2) y(j2,k0,m2) x(j1,m2) y(j1,k2,m2) y(j2,k1,m2) x(j2,m2) 
End
