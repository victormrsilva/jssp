\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,0) + x(0,2) >= 3
ord(0,2):  x(0,1) - x(0,2) >= 9
ord(1,1):  x(1,1) - x(1,2) >= 7
ord(1,2):  x(1,0) - x(1,1) >= 8
ord(2,1):  - x(2,0) + x(2,1) >= 7
ord(2,2):  - x(2,1) + x(2,2) >= 4
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 1
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 3
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -49 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  12 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 12
MCLin4(0,1,0):  12 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -49 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -49
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -16 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  45 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 45
MCLin9(1,0,0):  45 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -16 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -16
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 8
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 9
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -36 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  25 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 25
MCLin4(0,1,1):  25 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -36 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -36
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -42 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  19 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 19
MCLin9(1,0,1):  19 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -42 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -42
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 7
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 9
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -37 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  24 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 24
MCLin4(0,1,2):  24 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -37 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -37
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -40 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  21 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 21
MCLin9(1,0,2):  21 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -40 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -40
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 7
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 3
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -44 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  21 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 21
MCLin4(0,2,0):  21 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -44 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -44
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -31 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  34 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 34
MCLin9(2,0,0):  34 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -31 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -31
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 4
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 9
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -36 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  29 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 29
MCLin4(0,2,1):  29 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -36 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -36
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -42 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  23 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 23
MCLin9(2,0,1):  23 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -42 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -42
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 1
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 9
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -46 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  19 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 19
MCLin4(0,2,2):  19 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -46 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -46
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -29 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  36 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 36
MCLin9(2,0,2):  36 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -29 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -29
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 7
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 1
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -29 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  41 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 41
MCLin4(1,2,0):  41 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -29 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -29
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -49 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  21 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 21
MCLin9(2,1,0):  21 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -49 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -49
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 4
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 8
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -41 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  29 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 29
MCLin4(1,2,1):  29 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -41 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -41
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -41 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  29 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 29
MCLin9(2,1,1):  29 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -41 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -41
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 1
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 7
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -49 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  21 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 21
MCLin4(1,2,2):  21 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -49 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -49
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -29 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  41 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 41
MCLin9(2,1,2):  41 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -29 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -29
makespan(0):  C - x(0,1) >= 9
makespan(1):  C - x(1,0) >= 1
makespan(2):  C - x(2,2) >= 1
basic_cuts1(02,0):  3 x(0,0) + 7 x(2,0) >= 21
basic_cuts1(012,0):  3 x(0,0) + x(1,0) + 7 x(2,0) >= 31
cut_clique_1(02,0):  0.14286 x(0,0) + 0.33333 x(2,0) >= 1
cut_clique_1(012,0):  0.14286 x(0,0) + 0.33333 x(2,0) >= 1
basic_cuts_reverse1(01,1):  17 C -9 x(0,1) -8 x(1,1) >= 217
basic_cuts1(012,1):  9 x(0,1) + 8 x(1,1) + 4 x(2,1) >= 287
basic_cuts_reverse1(012,1):  21 C -9 x(0,1) -8 x(1,1) -4 x(2,1) >= 301
cut_clique_1(01,1):  0.06061 x(0,1) + 0.01299 x(1,1) >= 1
cut_clique_1(012,1):  0.03390 x(0,1) + 0.02373 x(1,1) + 0.01356 x(2,1) >= 1
basic_cuts1(01,2):  9 x(0,2) + 7 x(1,2) >= 63
basic_cuts_reverse1(01,2):  16 C -9 x(0,2) -7 x(1,2) >= 337
basic_cuts1(012,2):  9 x(0,2) + 7 x(1,2) + x(2,2) >= 79
cut_clique_1(01,2):  0.14286 x(0,2) + 0.04762 x(1,2) >= 1
cut_clique_1(012,2):  0.14286 x(0,2) + 0.04762 x(1,2) >= 1
cut_clique_6:  0.04348 x(0,2) + 0.02676 x(1,2) + 0.04348 x(2,2) >= 1
cut_clique_7:  0.12927 x(0,2) + 0.04160 x(1,2) + 0.00594 x(2,2) >= 1
Bounds
 0 <= C <= 49
 0 <= x(0,0) <= 28
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -49 <= v(0,1,0) <= 12
 -16 <= v(1,0,0) <= 45
 -490 <= u(0,1,0)
 -490 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -44 <= v(0,2,0) <= 21
 -31 <= v(2,0,0) <= 34
 -490 <= u(0,2,0)
 -490 <= u(2,0,0)
 15 <= x(1,0) <= 48
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -29 <= v(1,2,0) <= 41
 -49 <= v(2,1,0) <= 21
 -490 <= u(1,2,0)
 -490 <= u(2,1,0)
 0 <= x(2,0) <= 37
 12 <= x(0,1) <= 40
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -36 <= v(0,1,1) <= 25
 -42 <= v(1,0,1) <= 19
 -490 <= u(0,1,1)
 -490 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -36 <= v(0,2,1) <= 29
 -42 <= v(2,0,1) <= 23
 -490 <= u(0,2,1)
 -490 <= u(2,0,1)
 7 <= x(1,1) <= 40
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -41 <= v(1,2,1) <= 29
 -41 <= v(2,1,1) <= 29
 -490 <= u(1,2,1)
 -490 <= u(2,1,1)
 7 <= x(2,1) <= 44
 3 <= x(0,2) <= 31
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -37 <= v(0,1,2) <= 24
 -40 <= v(1,0,2) <= 21
 -490 <= u(0,1,2)
 -490 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -46 <= v(0,2,2) <= 19
 -29 <= v(2,0,2) <= 36
 -490 <= u(0,2,2)
 -490 <= u(2,0,2)
 0 <= x(1,2) <= 33
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -49 <= v(1,2,2) <= 21
 -29 <= v(2,1,2) <= 41
 -490 <= u(1,2,2)
 -490 <= u(2,1,2)
 11 <= x(2,2) <= 48
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) 
End
