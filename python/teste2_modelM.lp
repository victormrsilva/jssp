\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(j0,m0) + x(j0,m1) >= 10
ord(0,2):  - x(j0,m1) + x(j0,m2) >= 4
ord(1,1):  x(j1,m1) - x(j1,m2) >= 2
ord(1,2):  x(j1,m0) - x(j1,m1) >= 7
ord(2,1):  - x(j2,m1) + x(j2,m2) >= 10
ord(2,2):  x(j2,m0) - x(j2,m2) >= 1
phi(0,1,0):  x(j0,m0) + 54 y(j0,k1,m0) - x(j1,m0) >= 6
psy(0,1,0):  - x(j0,m0) -54 y(j0,k1,m0) + x(j1,m0) >= -44
phi(0,1,1):  x(j0,m1) + 54 y(j0,k1,m1) - x(j1,m1) >= 7
psy(0,1,1):  - x(j0,m1) -54 y(j0,k1,m1) + x(j1,m1) >= -50
phi(0,1,2):  x(j0,m2) + 54 y(j0,k1,m2) - x(j1,m2) >= 2
psy(0,1,2):  - x(j0,m2) -54 y(j0,k1,m2) + x(j1,m2) >= -50
phi(0,2,0):  x(j0,m0) + 54 y(j0,k2,m0) - x(j2,m0) >= 10
psy(0,2,0):  - x(j0,m0) -54 y(j0,k2,m0) + x(j2,m0) >= -44
phi(0,2,1):  x(j0,m1) + 54 y(j0,k2,m1) - x(j2,m1) >= 10
psy(0,2,1):  - x(j0,m1) -54 y(j0,k2,m1) + x(j2,m1) >= -50
phi(0,2,2):  x(j0,m2) + 54 y(j0,k2,m2) - x(j2,m2) >= 1
psy(0,2,2):  - x(j0,m2) -54 y(j0,k2,m2) + x(j2,m2) >= -50
phi(1,0,0):  - x(j0,m0) + 54 y(j1,k0,m0) + x(j1,m0) >= 10
psy(1,0,0):  x(j0,m0) -54 y(j1,k0,m0) - x(j1,m0) >= -48
phi(1,0,1):  - x(j0,m1) + 54 y(j1,k0,m1) + x(j1,m1) >= 4
psy(1,0,1):  x(j0,m1) -54 y(j1,k0,m1) - x(j1,m1) >= -47
phi(1,0,2):  - x(j0,m2) + 54 y(j1,k0,m2) + x(j1,m2) >= 4
psy(1,0,2):  x(j0,m2) -54 y(j1,k0,m2) - x(j1,m2) >= -52
phi(1,2,0):  x(j1,m0) + 54 y(j1,k2,m0) - x(j2,m0) >= 10
psy(1,2,0):  - x(j1,m0) -54 y(j1,k2,m0) + x(j2,m0) >= -48
phi(1,2,1):  x(j1,m1) + 54 y(j1,k2,m1) - x(j2,m1) >= 10
psy(1,2,1):  - x(j1,m1) -54 y(j1,k2,m1) + x(j2,m1) >= -47
phi(1,2,2):  x(j1,m2) + 54 y(j1,k2,m2) - x(j2,m2) >= 1
psy(1,2,2):  - x(j1,m2) -54 y(j1,k2,m2) + x(j2,m2) >= -52
phi(2,0,0):  - x(j0,m0) + 54 y(j2,k0,m0) + x(j2,m0) >= 10
psy(2,0,0):  x(j0,m0) -54 y(j2,k0,m0) - x(j2,m0) >= -44
phi(2,0,1):  - x(j0,m1) + 54 y(j2,k0,m1) + x(j2,m1) >= 4
psy(2,0,1):  x(j0,m1) -54 y(j2,k0,m1) - x(j2,m1) >= -44
phi(2,0,2):  - x(j0,m2) + 54 y(j2,k0,m2) + x(j2,m2) >= 4
psy(2,0,2):  x(j0,m2) -54 y(j2,k0,m2) - x(j2,m2) >= -53
phi(2,1,0):  - x(j1,m0) + 54 y(j2,k1,m0) + x(j2,m0) >= 6
psy(2,1,0):  x(j1,m0) -54 y(j2,k1,m0) - x(j2,m0) >= -44
phi(2,1,1):  - x(j1,m1) + 54 y(j2,k1,m1) + x(j2,m1) >= 7
psy(2,1,1):  x(j1,m1) -54 y(j2,k1,m1) - x(j2,m1) >= -44
phi(2,1,2):  - x(j1,m2) + 54 y(j2,k1,m2) + x(j2,m2) >= 2
psy(2,1,2):  x(j1,m2) -54 y(j2,k1,m2) - x(j2,m2) >= -53
makespan(0):  C - x(j0,m2) >= 4
makespan(1):  C - x(j1,m0) >= 6
makespan(2):  C - x(j2,m0) >= 10
Bounds
 0 <= C <= 54
 0 <= x(j0,m0) <= 36
 0 <= y(j0,k1,m0) <= 1
 0 <= y(j1,k0,m0) <= 1
 0 <= y(j0,k2,m0) <= 1
 0 <= y(j2,k0,m0) <= 1
 9 <= x(j1,m0) <= 48
 0 <= y(j1,k2,m0) <= 1
 0 <= y(j2,k1,m0) <= 1
 11 <= x(j2,m0) <= 44
 10 <= x(j0,m1) <= 46
 0 <= y(j0,k1,m1) <= 1
 0 <= y(j1,k0,m1) <= 1
 0 <= y(j0,k2,m1) <= 1
 0 <= y(j2,k0,m1) <= 1
 2 <= x(j1,m1) <= 41
 0 <= y(j1,k2,m1) <= 1
 0 <= y(j2,k1,m1) <= 1
 0 <= x(j2,m1) <= 33
 14 <= x(j0,m2) <= 50
 0 <= y(j0,k1,m2) <= 1
 0 <= y(j1,k0,m2) <= 1
 0 <= y(j0,k2,m2) <= 1
 0 <= y(j2,k0,m2) <= 1
 0 <= x(j1,m2) <= 39
 0 <= y(j1,k2,m2) <= 1
 0 <= y(j2,k1,m2) <= 1
 10 <= x(j2,m2) <= 43
Integers
C x(j0,m0) y(j0,k1,m0) y(j1,k0,m0) y(j0,k2,m0) y(j2,k0,m0) x(j1,m0) y(j1,k2,m0) y(j2,k1,m0) x(j2,m0)
x(j0,m1) y(j0,k1,m1) y(j1,k0,m1) y(j0,k2,m1) y(j2,k0,m1) x(j1,m1) y(j1,k2,m1) y(j2,k1,m1) x(j2,m1) x(j0,m2)
y(j0,k1,m2) y(j1,k0,m2) y(j0,k2,m2) y(j2,k0,m2) x(j1,m2) y(j1,k2,m2) y(j2,k1,m2) x(j2,m2)
End
