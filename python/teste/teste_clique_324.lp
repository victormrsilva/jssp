\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,0) + x(0,1) >= 10
ord(0,2):  - x(0,1) + x(0,2) >= 3
ord(1,1):  x(1,1) - x(1,2) >= 1
ord(1,2):  x(1,0) - x(1,1) >= 4
ord(2,1):  x(2,0) - x(2,1) >= 9
ord(2,2):  - x(2,0) + x(2,2) >= 7
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 2
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 10
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -51 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  22 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 22
MCLin4(0,1,0):  22 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -51 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -51
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -34 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  39 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 39
MCLin9(1,0,0):  39 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -34 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -34
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 4
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 3
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -39 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  34 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 34
MCLin4(0,1,1):  34 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -39 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -39
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -41 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  32 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 32
MCLin9(1,0,1):  32 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -41 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -41
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 1
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 9
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -32 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  41 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 41
MCLin4(0,1,2):  41 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -32 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -32
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -51 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  22 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 22
MCLin9(1,0,2):  22 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -51 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -51
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 7
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 10
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -45 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  13 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 13
MCLin4(0,2,0):  13 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -45 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -45
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -30 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  28 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 28
MCLin9(2,0,0):  28 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -30 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -30
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 9
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 3
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -28 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  30 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 30
MCLin4(0,2,1):  30 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -28 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -28
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -42 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  16 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 16
MCLin9(2,0,1):  16 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -42 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -42
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 6
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 9
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -38 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  20 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 20
MCLin4(0,2,2):  20 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -38 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -38
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -35 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  23 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 23
MCLin9(2,0,2):  23 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -35 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -35
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 7
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 2
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -40 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  33 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 33
MCLin4(1,2,0):  33 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -40 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -40
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -42 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  31 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 31
MCLin9(2,1,0):  31 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -42 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -42
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 9
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 4
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -37 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  36 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 36
MCLin4(1,2,1):  36 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -37 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -37
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -49 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  24 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 24
MCLin9(2,1,1):  24 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -49 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -49
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 6
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 1
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -51 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  22 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 22
MCLin4(1,2,2):  22 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -51 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -51
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -29 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  44 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 44
MCLin9(2,1,2):  44 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -29 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -29
makespan(0):  C - x(0,2) >= 9
makespan(1):  C - x(1,0) >= 2
makespan(2):  C - x(2,2) >= 6
cut_clique_0(012,0):  0.01863 x(0,0) + 0.01923 x(1,0) + 0.06731 x(2,0) >= 1
cut_clique_0(012,1):  0.11111 x(1,1) + 0.17778 x(2,1) >= 1
cut_clique_0(012,2):  0.02174 x(0,2) + 0.03261 x(2,2) >= 1
cut_clique_1(012,0):  0.07143 x(0,0) + 0.10000 x(1,0) >= 1
cut_clique_1(012,1):  0.05868 x(0,1) + 0.03178 x(1,1) + 0.02934 x(2,1) >= 1
cut_clique_0:  0.00290 x(0,0) + 0.01126 x(1,0) + 0.00901 x(2,0) + 0.04639 x(0,1) + 0.01959 x(1,1) + 0.02207 x(2,1) >= 1
cut_clique_1:  0.01418 x(0,0) + 0.01945 x(1,0) + 0.05186 x(2,0) + 0.01297 x(0,1) + 0.00162 x(1,1) >= 1
Bounds
 0 <= C <= 51
 0 <= x(0,0) <= 29
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -51 <= v(0,1,0) <= 22
 -34 <= v(1,0,0) <= 39
 -510 <= u(0,1,0)
 -510 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -45 <= v(0,2,0) <= 13
 -30 <= v(2,0,0) <= 28
 -510 <= u(0,2,0)
 -510 <= u(2,0,0)
 5 <= x(1,0) <= 49
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -40 <= v(1,2,0) <= 33
 -42 <= v(2,1,0) <= 31
 -510 <= u(1,2,0)
 -510 <= u(2,1,0)
 9 <= x(2,0) <= 38
 10 <= x(0,1) <= 39
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -39 <= v(0,1,1) <= 34
 -41 <= v(1,0,1) <= 32
 -510 <= u(0,1,1)
 -510 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -28 <= v(0,2,1) <= 30
 -42 <= v(2,0,1) <= 16
 -510 <= u(0,2,1)
 -510 <= u(2,0,1)
 1 <= x(1,1) <= 45
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -37 <= v(1,2,1) <= 36
 -49 <= v(2,1,1) <= 24
 -510 <= u(1,2,1)
 -510 <= u(2,1,1)
 0 <= x(2,1) <= 29
 13 <= x(0,2) <= 42
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -32 <= v(0,1,2) <= 41
 -51 <= v(1,0,2) <= 22
 -510 <= u(0,1,2)
 -510 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -38 <= v(0,2,2) <= 20
 -35 <= v(2,0,2) <= 23
 -510 <= u(0,2,2)
 -510 <= u(2,0,2)
 0 <= x(1,2) <= 44
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -51 <= v(1,2,2) <= 22
 -29 <= v(2,1,2) <= 44
 -510 <= u(1,2,2)
 -510 <= u(2,1,2)
 16 <= x(2,2) <= 45
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) 
End
