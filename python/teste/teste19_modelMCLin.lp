\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  x(0,0) - x(0,1) >= 4
ord(0,2):  - x(0,0) + x(0,2) >= 6
ord(1,1):  - x(1,0) + x(1,2) >= 5
ord(1,2):  x(1,1) - x(1,2) >= 1
ord(2,1):  - x(2,0) + x(2,2) >= 9
ord(2,2):  x(2,1) - x(2,2) >= 10
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 5
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 6
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -40 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  40 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 40
MCLin4(0,1,0):  40 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -40 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -40
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -51 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  29 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 29
MCLin9(1,0,0):  29 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -51 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -51
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 8
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 4
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -53 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  27 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 27
MCLin4(0,1,1):  27 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -53 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -53
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -39 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  41 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 41
MCLin9(1,0,1):  41 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -39 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -39
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 1
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 2
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -35 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  45 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 45
MCLin4(0,1,2):  45 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -35 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -35
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -48 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  32 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 32
MCLin9(1,0,2):  32 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -48 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -48
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 9
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 6
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -31 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  36 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 36
MCLin4(0,2,0):  36 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -31 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -31
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -51 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  16 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 16
MCLin9(2,0,0):  16 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -51 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -51
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 8
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 4
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -53 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  14 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 14
MCLin4(0,2,1):  14 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -53 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -53
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -26 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  41 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 41
MCLin9(2,0,1):  41 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -26 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -26
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 10
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 2
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -35 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  32 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 32
MCLin4(0,2,2):  32 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -35 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -35
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -44 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  23 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 23
MCLin9(2,0,2):  23 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -44 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -44
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 9
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 5
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -35 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  30 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 30
MCLin4(1,2,0):  30 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -35 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -35
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -44 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  21 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 21
MCLin9(2,1,0):  21 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -44 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -44
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 8
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 8
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -47 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  18 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 18
MCLin4(1,2,1):  18 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -47 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -47
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -34 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  31 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 31
MCLin9(2,1,1):  31 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -34 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -34
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 10
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 1
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -40 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  25 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 25
MCLin4(1,2,2):  25 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -40 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -40
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -36 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  29 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 29
MCLin9(2,1,2):  29 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -36 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -36
makespan(0):  C - x(0,2) >= 2
makespan(1):  C - x(1,1) >= 8
makespan(2):  C - x(2,1) >= 8
Bounds
 0 <= C <= 53
 4 <= x(0,0) <= 45
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -40 <= v(0,1,0) <= 40
 -51 <= v(1,0,0) <= 29
 -530 <= u(0,1,0)
 -530 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -31 <= v(0,2,0) <= 36
 -51 <= v(2,0,0) <= 16
 -530 <= u(0,2,0)
 -530 <= u(2,0,0)
 0 <= x(1,0) <= 39
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -35 <= v(1,2,0) <= 30
 -44 <= v(2,1,0) <= 21
 -530 <= u(1,2,0)
 -530 <= u(2,1,0)
 0 <= x(2,0) <= 26
 0 <= x(0,1) <= 41
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -53 <= v(0,1,1) <= 27
 -39 <= v(1,0,1) <= 41
 -530 <= u(0,1,1)
 -530 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -53 <= v(0,2,1) <= 14
 -26 <= v(2,0,1) <= 41
 -530 <= u(0,2,1)
 -530 <= u(2,0,1)
 6 <= x(1,1) <= 45
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -47 <= v(1,2,1) <= 18
 -34 <= v(2,1,1) <= 31
 -530 <= u(1,2,1)
 -530 <= u(2,1,1)
 19 <= x(2,1) <= 45
 10 <= x(0,2) <= 51
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -35 <= v(0,1,2) <= 45
 -48 <= v(1,0,2) <= 32
 -530 <= u(0,1,2)
 -530 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -35 <= v(0,2,2) <= 32
 -44 <= v(2,0,2) <= 23
 -530 <= u(0,2,2)
 -530 <= u(2,0,2)
 5 <= x(1,2) <= 44
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -40 <= v(1,2,2) <= 25
 -36 <= v(2,1,2) <= 29
 -530 <= u(1,2,2)
 -530 <= u(2,1,2)
 9 <= x(2,2) <= 35
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) 
End
