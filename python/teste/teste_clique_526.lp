\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  x(0,0) - x(0,1) >= 4
ord(0,2):  - x(0,0) + x(0,2) >= 2
ord(1,1):  x(1,0) - x(1,1) >= 2
ord(1,2):  - x(1,0) + x(1,2) >= 2
ord(2,1):  x(2,1) - x(2,2) >= 1
ord(2,2):  x(2,0) - x(2,1) >= 4
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
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 2
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -22 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  20 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 20
MCLin4(0,1,0):  20 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -22 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -22
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -24 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  18 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 18
MCLin9(1,0,0):  18 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -24 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -24
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 2
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 4
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -24 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  18 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 18
MCLin4(0,1,1):  18 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -24 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -24
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -24 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  18 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 18
MCLin9(1,0,1):  18 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -24 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -24
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 4
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 4
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -24 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  18 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 18
MCLin4(0,1,2):  18 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -24 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -24
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -26 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  16 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 16
MCLin9(1,0,2):  16 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -26 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -26
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 7
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 2
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -26 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  12 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 12
MCLin4(0,2,0):  12 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -26 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -26
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -21 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  17 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 17
MCLin9(2,0,0):  17 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -21 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -21
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 4
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 4
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -23 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  15 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 15
MCLin4(0,2,1):  15 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -23 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -23
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -23 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  15 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 15
MCLin9(2,0,1):  15 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -23 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -23
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 1
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 4
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -13 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  25 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 25
MCLin4(0,2,2):  25 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -13 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -13
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -30 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  8 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 8
MCLin9(2,0,2):  8 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -30 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -30
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 7
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 2
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -28 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  12 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 12
MCLin4(1,2,0):  12 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -28 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -28
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -21 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  19 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 19
MCLin9(2,1,0):  19 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -21 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -21
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 4
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 2
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -23 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  17 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 17
MCLin4(1,2,1):  17 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -23 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -23
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -23 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  17 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 17
MCLin9(2,1,1):  17 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -23 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -23
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 1
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 4
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -15 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  25 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 25
MCLin4(1,2,2):  25 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -15 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -15
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -30 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  10 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 10
MCLin9(2,1,2):  10 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -30 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -30
makespan(0):  C - x(0,2) >= 4
makespan(1):  C - x(1,2) >= 4
makespan(2):  C - x(2,0) >= 7
cut_clique_0(012,0):  0.01923 x(0,0) + 0.15385 x(2,0) >= 1
cut_clique_0(012,1):  0.50000 x(0,1) + 0.25000 x(1,1) >= 1
cut_clique_0(012,2):  0.10714 x(0,2) + 0.03571 x(1,2) >= 1
cut_clique_1(012,1):  0.12500 x(0,1) + 0.06250 x(1,1) + 0.12500 x(2,1) >= 1
cut_clique_2(012,1):  0.10484 x(0,1) + 0.03226 x(1,1) + 0.18548 x(2,1) >= 1
cut_clique_0:  0.00973 x(0,0) + 0.01946 x(1,0) + 0.02433 x(2,0) + 0.07360 x(0,1) + 0.14051 x(2,1) >= 1
cut_clique_1:  0.01256 x(0,0) + 0.02010 x(1,0) + 0.02513 x(2,0) + 0.06910 x(0,1) + 0.13442 x(2,1) >= 1
Bounds
 0 <= C <= 30
 4 <= x(0,0) <= 24
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -22 <= v(0,1,0) <= 20
 -24 <= v(1,0,0) <= 18
 -300 <= u(0,1,0)
 -300 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -26 <= v(0,2,0) <= 12
 -21 <= v(2,0,0) <= 17
 -300 <= u(0,2,0)
 -300 <= u(2,0,0)
 2 <= x(1,0) <= 24
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -28 <= v(1,2,0) <= 12
 -21 <= v(2,1,0) <= 19
 -300 <= u(1,2,0)
 -300 <= u(2,1,0)
 5 <= x(2,0) <= 23
 0 <= x(0,1) <= 20
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -24 <= v(0,1,1) <= 18
 -24 <= v(1,0,1) <= 18
 -300 <= u(0,1,1)
 -300 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -23 <= v(0,2,1) <= 15
 -23 <= v(2,0,1) <= 15
 -300 <= u(0,2,1)
 -300 <= u(2,0,1)
 0 <= x(1,1) <= 22
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -23 <= v(1,2,1) <= 17
 -23 <= v(2,1,1) <= 17
 -300 <= u(1,2,1)
 -300 <= u(2,1,1)
 1 <= x(2,1) <= 19
 6 <= x(0,2) <= 26
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -24 <= v(0,1,2) <= 18
 -26 <= v(1,0,2) <= 16
 -300 <= u(0,1,2)
 -300 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -13 <= v(0,2,2) <= 25
 -30 <= v(2,0,2) <= 8
 -300 <= u(0,2,2)
 -300 <= u(2,0,2)
 4 <= x(1,2) <= 26
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -15 <= v(1,2,2) <= 25
 -30 <= v(2,1,2) <= 10
 -300 <= u(1,2,2)
 -300 <= u(2,1,2)
 0 <= x(2,2) <= 18
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) 
u(2,1,0) x(2,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) 
y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) 
u(1,2,1) u(2,1,1) x(2,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) 
y(0,2,2) y(2,0,2) v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) x(2,2) 
End
