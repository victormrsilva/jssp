\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,3) >= 10
ord(0,2):  x(0,0) - x(0,3) >= 3
ord(0,3):  - x(0,0) + x(0,2) >= 2
ord(1,1):  - x(1,1) + x(1,3) >= 2
ord(1,2):  x(1,0) - x(1,3) >= 10
ord(1,3):  - x(1,0) + x(1,2) >= 9
ord(2,1):  x(2,0) - x(2,2) >= 7
ord(2,2):  - x(2,0) + x(2,3) >= 3
ord(2,3):  x(2,1) - x(2,3) >= 5
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,1,3):  y(0,1,3) + y(1,0,3) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(0,2,3):  y(0,2,3) + y(2,0,3) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
triangle2(1,2,3):  y(1,2,3) + y(2,1,3) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 9
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 2
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -51 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  44 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 44
MCLin4(0,1,0):  44 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -51 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -51
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -55 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  40 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 40
MCLin9(1,0,0):  40 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -55 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -55
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 2
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 10
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -45 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  50 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 50
MCLin4(0,1,1):  50 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -45 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -45
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -62 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  33 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 33
MCLin9(1,0,1):  33 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -62 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -62
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 10
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 7
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -59 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  36 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 36
MCLin4(0,1,2):  36 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -59 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -59
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -53 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  42 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 42
MCLin9(1,0,2):  42 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -53 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -53
MCLinV(0,1,3):  x(0,3) - v(0,1,3) - x(1,3) = 10
MCLinV(1,0,3):  - x(0,3) - v(1,0,3) + x(1,3) = 3
MCLin1(0,1,3):  u(0,1,3) >= -0
MCLin2(0,1,3):  -45 y(1,0,3) - u(0,1,3) <= -0
MCLin3(0,1,3):  50 y(1,0,3) + v(0,1,3) - u(0,1,3) <= 50
MCLin4(0,1,3):  50 y(1,0,3) - u(0,1,3) >= -0
MCLin5(0,1,3):  -45 y(1,0,3) + v(0,1,3) - u(0,1,3) >= -45
MCLin6(1,0,3):  u(1,0,3) >= -0
MCLin7(1,0,3):  -63 y(0,1,3) - u(1,0,3) <= -0
MCLin8(1,0,3):  32 y(0,1,3) + v(1,0,3) - u(1,0,3) <= 32
MCLin9(1,0,3):  32 y(0,1,3) - u(1,0,3) >= -0
MCLin10(1,0,3):  -63 y(0,1,3) + v(1,0,3) - u(1,0,3) >= -63
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 3
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 2
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -50 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  55 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 55
MCLin4(0,2,0):  55 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -50 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -50
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -60 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  45 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 45
MCLin9(2,0,0):  45 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -60 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -60
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 6
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 10
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -74 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  31 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 31
MCLin4(0,2,1):  31 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -74 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -74
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -47 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  58 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 58
MCLin9(2,0,1):  58 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -47 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -47
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 7
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 7
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -45 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  60 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 60
MCLin4(0,2,2):  60 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -45 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -45
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -74 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  31 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 31
MCLin9(2,0,2):  31 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -74 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -74
MCLinV(0,2,3):  x(0,3) - v(0,2,3) - x(2,3) = 5
MCLinV(2,0,3):  - x(0,3) - v(2,0,3) + x(2,3) = 3
MCLin1(0,2,3):  u(0,2,3) >= -0
MCLin2(0,2,3):  -58 y(2,0,3) - u(0,2,3) <= -0
MCLin3(0,2,3):  47 y(2,0,3) + v(0,2,3) - u(0,2,3) <= 47
MCLin4(0,2,3):  47 y(2,0,3) - u(0,2,3) >= -0
MCLin5(0,2,3):  -58 y(2,0,3) + v(0,2,3) - u(0,2,3) >= -58
MCLin6(2,0,3):  u(2,0,3) >= -0
MCLin7(2,0,3):  -55 y(0,2,3) - u(2,0,3) <= -0
MCLin8(2,0,3):  50 y(0,2,3) + v(2,0,3) - u(2,0,3) <= 50
MCLin9(2,0,3):  50 y(0,2,3) - u(2,0,3) >= -0
MCLin10(2,0,3):  -55 y(0,2,3) + v(2,0,3) - u(2,0,3) >= -55
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 3
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 9
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -51 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  45 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 45
MCLin4(1,2,0):  45 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -51 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -51
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -57 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  39 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 39
MCLin9(2,1,0):  39 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -57 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -57
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 6
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 2
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -74 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  22 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 22
MCLin4(1,2,1):  22 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -74 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -74
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -30 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  66 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 66
MCLin9(2,1,1):  66 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -30 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -30
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 7
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 10
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -39 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  57 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 57
MCLin4(1,2,2):  57 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -39 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -39
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -74 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  22 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 22
MCLin9(2,1,2):  22 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -74 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -74
MCLinV(1,2,3):  x(1,3) - v(1,2,3) - x(2,3) = 5
MCLinV(2,1,3):  - x(1,3) - v(2,1,3) + x(2,3) = 10
MCLin1(1,2,3):  u(1,2,3) >= -0
MCLin2(1,2,3):  -66 y(2,1,3) - u(1,2,3) <= -0
MCLin3(1,2,3):  30 y(2,1,3) + v(1,2,3) - u(1,2,3) <= 30
MCLin4(1,2,3):  30 y(2,1,3) - u(1,2,3) >= -0
MCLin5(1,2,3):  -66 y(2,1,3) + v(1,2,3) - u(1,2,3) >= -66
MCLin6(2,1,3):  u(2,1,3) >= -0
MCLin7(2,1,3):  -45 y(1,2,3) - u(2,1,3) <= -0
MCLin8(2,1,3):  51 y(1,2,3) + v(2,1,3) - u(2,1,3) <= 51
MCLin9(2,1,3):  51 y(1,2,3) - u(2,1,3) >= -0
MCLin10(2,1,3):  -45 y(1,2,3) + v(2,1,3) - u(2,1,3) >= -45
makespan(0):  C - x(0,2) >= 7
makespan(1):  C - x(1,2) >= 10
makespan(2):  C - x(2,1) >= 6
cut_clique_0(012,0):  0.01887 x(0,0) + 0.05031 x(1,0) >= 1
cut_clique_0(012,2):  0.00272 x(0,2) + 0.04360 x(1,2) >= 1
cut_clique_0(012,3):  0.02655 x(0,3) + 0.00885 x(1,3) + 0.04425 x(2,3) >= 1
cut_clique_1(012,3):  0.01250 x(1,3) + 0.08125 x(2,3) >= 1
cut_clique_2(012,3):  0.00302 x(0,3) + 0.01005 x(1,3) + 0.07739 x(2,3) >= 1
cut_clique_0:  0.00566 x(0,2) + 0.02641 x(1,2) + 0.00668 x(1,3) + 0.02138 x(2,3) >= 1
cut_clique_1:  0.04072 x(1,0) + 0.00838 x(0,2) + 0.01198 x(1,2) >= 1
cut_clique_2:  0.01148 x(1,0) + 0.00827 x(0,2) + 0.01181 x(1,2) + 0.00681 x(1,3) + 0.02869 x(2,3) >= 1
cut_clique_3:  0.01887 x(0,0) + 0.02673 x(1,0) + 0.02358 x(2,3) >= 1
cut_clique_4:  0.02204 x(1,0) + 0.00454 x(0,2) + 0.00648 x(1,2) + 0.01240 x(0,3) + 0.02066 x(2,3) >= 1
cut_clique_5:  0.02545 x(1,0) + 0.00838 x(0,2) + 0.01198 x(1,2) + 0.01527 x(2,3) >= 1
Bounds
 0 <= C <= 74
 13 <= x(0,0) <= 65
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -51 <= v(0,1,0) <= 44
 -55 <= v(1,0,0) <= 40
 -740 <= u(0,1,0)
 -740 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -50 <= v(0,2,0) <= 55
 -60 <= v(2,0,0) <= 45
 -740 <= u(0,2,0)
 -740 <= u(2,0,0)
 12 <= x(1,0) <= 55
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -51 <= v(1,2,0) <= 45
 -57 <= v(2,1,0) <= 39
 -740 <= u(1,2,0)
 -740 <= u(2,1,0)
 7 <= x(2,0) <= 60
 0 <= x(0,1) <= 52
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -45 <= v(0,1,1) <= 50
 -62 <= v(1,0,1) <= 33
 -740 <= u(0,1,1)
 -740 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -74 <= v(0,2,1) <= 31
 -47 <= v(2,0,1) <= 58
 -740 <= u(0,2,1)
 -740 <= u(2,0,1)
 0 <= x(1,1) <= 43
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -74 <= v(1,2,1) <= 22
 -30 <= v(2,1,1) <= 66
 -740 <= u(1,2,1)
 -740 <= u(2,1,1)
 15 <= x(2,1) <= 68
 15 <= x(0,2) <= 67
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -59 <= v(0,1,2) <= 36
 -53 <= v(1,0,2) <= 42
 -740 <= u(0,1,2)
 -740 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -45 <= v(0,2,2) <= 60
 -74 <= v(2,0,2) <= 31
 -740 <= u(0,2,2)
 -740 <= u(2,0,2)
 21 <= x(1,2) <= 64
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -39 <= v(1,2,2) <= 57
 -74 <= v(2,1,2) <= 22
 -740 <= u(1,2,2)
 -740 <= u(2,1,2)
 0 <= x(2,2) <= 53
 10 <= x(0,3) <= 62
 0 <= y(0,1,3) <= 1
 0 <= y(1,0,3) <= 1
 -45 <= v(0,1,3) <= 50
 -63 <= v(1,0,3) <= 32
 -740 <= u(0,1,3)
 -740 <= u(1,0,3)
 0 <= y(0,2,3) <= 1
 0 <= y(2,0,3) <= 1
 -58 <= v(0,2,3) <= 47
 -55 <= v(2,0,3) <= 50
 -740 <= u(0,2,3)
 -740 <= u(2,0,3)
 2 <= x(1,3) <= 45
 0 <= y(1,2,3) <= 1
 0 <= y(2,1,3) <= 1
 -66 <= v(1,2,3) <= 30
 -45 <= v(2,1,3) <= 51
 -740 <= u(1,2,3)
 -740 <= u(2,1,3)
 10 <= x(2,3) <= 63
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) x(0,3) y(0,1,3) y(1,0,3) v(0,1,3) v(1,0,3) u(0,1,3) 
u(1,0,3) y(0,2,3) y(2,0,3) v(0,2,3) v(2,0,3) u(0,2,3) u(2,0,3) x(1,3) y(1,2,3) y(2,1,3) 
v(1,2,3) v(2,1,3) u(1,2,3) u(2,1,3) x(2,3) 
End
