\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,0) + x(0,2) >= 5
ord(0,2):  x(0,1) - x(0,2) >= 6
ord(1,1):  x(1,0) - x(1,1) >= 9
ord(1,2):  - x(1,0) + x(1,2) >= 1
ord(2,1):  - x(2,1) + x(2,2) >= 2
ord(2,2):  x(2,0) - x(2,2) >= 2
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
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 5
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -34 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  20 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 20
MCLin4(0,1,0):  20 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -34 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -34
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -26 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  28 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 28
MCLin9(1,0,0):  28 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -26 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -26
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 9
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 1
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -22 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  32 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 32
MCLin4(0,1,1):  32 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -22 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -22
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -42 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  12 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 12
MCLin9(1,0,1):  12 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -42 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -42
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 8
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 6
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -37 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  17 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 17
MCLin4(0,1,2):  17 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -37 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -37
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -31 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  23 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 23
MCLin9(1,0,2):  23 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -31 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -31
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 8
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 5
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -42 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  18 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 18
MCLin4(0,2,0):  18 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -42 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -42
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -31 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  29 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 29
MCLin9(2,0,0):  29 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -31 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -31
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 1
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -21 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  39 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 39
MCLin4(0,2,1):  39 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -21 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -21
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -42 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  18 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 18
MCLin9(2,0,1):  18 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -42 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -42
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 2
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 6
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -29 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  31 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 31
MCLin4(0,2,2):  31 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -29 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -29
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -39 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  21 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 21
MCLin9(2,0,2):  21 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -39 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -39
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 8
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 1
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -33 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  21 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 21
MCLin4(1,2,0):  21 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -33 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -33
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -30 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  24 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 24
MCLin9(2,1,0):  24 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -30 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -30
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 9
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -32 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  22 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 22
MCLin4(1,2,1):  22 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -32 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -32
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -33 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  21 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 21
MCLin9(2,1,1):  21 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -33 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -33
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 2
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 8
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -24 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  30 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 30
MCLin4(1,2,2):  30 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -24 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -24
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -40 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  14 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 14
MCLin9(2,1,2):  14 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -40 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -40
makespan(0):  C - x(0,1) >= 1
makespan(1):  C - x(1,2) >= 8
makespan(2):  C - x(2,0) >= 8
Bounds
 0 <= C <= 42
 0 <= x(0,0) <= 30
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -34 <= v(0,1,0) <= 20
 -26 <= v(1,0,0) <= 28
 -420 <= u(0,1,0)
 -420 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -42 <= v(0,2,0) <= 18
 -31 <= v(2,0,0) <= 29
 -420 <= u(0,2,0)
 -420 <= u(2,0,0)
 9 <= x(1,0) <= 33
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -33 <= v(1,2,0) <= 21
 -30 <= v(2,1,0) <= 24
 -420 <= u(1,2,0)
 -420 <= u(2,1,0)
 4 <= x(2,0) <= 34
 11 <= x(0,1) <= 41
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -22 <= v(0,1,1) <= 32
 -42 <= v(1,0,1) <= 12
 -420 <= u(0,1,1)
 -420 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -21 <= v(0,2,1) <= 39
 -42 <= v(2,0,1) <= 18
 -420 <= u(0,2,1)
 -420 <= u(2,0,1)
 0 <= x(1,1) <= 24
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -32 <= v(1,2,1) <= 22
 -33 <= v(2,1,1) <= 21
 -420 <= u(1,2,1)
 -420 <= u(2,1,1)
 0 <= x(2,1) <= 30
 5 <= x(0,2) <= 35
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -37 <= v(0,1,2) <= 17
 -31 <= v(1,0,2) <= 23
 -420 <= u(0,1,2)
 -420 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -29 <= v(0,2,2) <= 31
 -39 <= v(2,0,2) <= 21
 -420 <= u(0,2,2)
 -420 <= u(2,0,2)
 10 <= x(1,2) <= 34
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -24 <= v(1,2,2) <= 30
 -40 <= v(2,1,2) <= 14
 -420 <= u(1,2,2)
 -420 <= u(2,1,2)
 2 <= x(2,2) <= 32
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) 
End
