\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,2) >= 5
ord(0,2):  - x(0,2) + x(0,3) >= 2
ord(0,3):  x(0,0) - x(0,3) >= 10
ord(1,1):  - x(1,1) + x(1,3) >= 8
ord(1,2):  x(1,2) - x(1,3) >= 10
ord(1,3):  x(1,0) - x(1,2) >= 1
ord(2,1):  - x(2,1) + x(2,2) >= 2
ord(2,2):  - x(2,2) + x(2,3) >= 3
ord(2,3):  x(2,0) - x(2,3) >= 1
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
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 4
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -45 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  30 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 30
MCLin4(0,1,0):  30 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -45 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -45
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -43 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  32 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 32
MCLin9(1,0,0):  32 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -43 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -43
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 8
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 5
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -42 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  33 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 33
MCLin4(0,1,1):  33 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -42 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -42
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -46 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  29 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 29
MCLin9(1,0,1):  29 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -46 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -46
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 1
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 2
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -48 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  27 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 27
MCLin4(0,1,2):  27 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -48 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -48
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -30 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  45 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 45
MCLin9(1,0,2):  45 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -30 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -30
MCLinV(0,1,3):  x(0,3) - v(0,1,3) - x(1,3) = 10
MCLinV(1,0,3):  - x(0,3) - v(1,0,3) + x(1,3) = 10
MCLin1(0,1,3):  u(0,1,3) >= -0
MCLin2(0,1,3):  -45 y(1,0,3) - u(0,1,3) <= -0
MCLin3(0,1,3):  30 y(1,0,3) + v(0,1,3) - u(0,1,3) <= 30
MCLin4(0,1,3):  30 y(1,0,3) - u(0,1,3) >= -0
MCLin5(0,1,3):  -45 y(1,0,3) + v(0,1,3) - u(0,1,3) >= -45
MCLin6(1,0,3):  u(1,0,3) >= -0
MCLin7(1,0,3):  -50 y(0,1,3) - u(1,0,3) <= -0
MCLin8(1,0,3):  25 y(0,1,3) + v(1,0,3) - u(1,0,3) <= 25
MCLin9(1,0,3):  25 y(0,1,3) - u(1,0,3) >= -0
MCLin10(1,0,3):  -50 y(0,1,3) + v(1,0,3) - u(1,0,3) >= -50
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 7
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 4
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -45 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  45 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 45
MCLin4(0,2,0):  45 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -45 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -45
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -56 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  34 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 34
MCLin9(2,0,0):  34 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -56 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -56
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 5
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -51 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  39 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 39
MCLin4(0,2,1):  39 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -51 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -51
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -46 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  44 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 44
MCLin9(2,0,1):  44 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -46 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -46
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 3
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 2
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -49 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  41 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 41
MCLin4(0,2,2):  41 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -49 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -49
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -46 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  44 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 44
MCLin9(2,0,2):  44 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -46 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -46
MCLinV(0,2,3):  x(0,3) - v(0,2,3) - x(2,3) = 1
MCLinV(2,0,3):  - x(0,3) - v(2,0,3) + x(2,3) = 10
MCLin1(0,2,3):  u(0,2,3) >= -0
MCLin2(0,2,3):  -48 y(2,0,3) - u(0,2,3) <= -0
MCLin3(0,2,3):  42 y(2,0,3) + v(0,2,3) - u(0,2,3) <= 42
MCLin4(0,2,3):  42 y(2,0,3) - u(0,2,3) >= -0
MCLin5(0,2,3):  -48 y(2,0,3) + v(0,2,3) - u(0,2,3) >= -48
MCLin6(2,0,3):  u(2,0,3) >= -0
MCLin7(2,0,3):  -53 y(0,2,3) - u(2,0,3) <= -0
MCLin8(2,0,3):  37 y(0,2,3) + v(2,0,3) - u(2,0,3) <= 37
MCLin9(2,0,3):  37 y(0,2,3) - u(2,0,3) >= -0
MCLin10(2,0,3):  -53 y(0,2,3) + v(2,0,3) - u(2,0,3) >= -53
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 7
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 9
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -43 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  40 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 40
MCLin4(1,2,0):  40 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -43 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -43
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -56 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  27 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 27
MCLin9(2,1,0):  27 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -56 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -56
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 8
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -51 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  32 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 32
MCLin4(1,2,1):  32 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -51 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -51
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -42 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  41 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 41
MCLin9(2,1,1):  41 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -42 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -42
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 3
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 1
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -36 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  47 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 47
MCLin4(1,2,2):  47 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -36 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -36
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -51 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  32 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 32
MCLin9(2,1,2):  32 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -51 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -51
MCLinV(1,2,3):  x(1,3) - v(1,2,3) - x(2,3) = 1
MCLinV(2,1,3):  - x(1,3) - v(2,1,3) + x(2,3) = 10
MCLin1(1,2,3):  u(1,2,3) >= -0
MCLin2(1,2,3):  -47 y(2,1,3) - u(1,2,3) <= -0
MCLin3(1,2,3):  36 y(2,1,3) + v(1,2,3) - u(1,2,3) <= 36
MCLin4(1,2,3):  36 y(2,1,3) - u(1,2,3) >= -0
MCLin5(1,2,3):  -47 y(2,1,3) + v(1,2,3) - u(1,2,3) >= -47
MCLin6(2,1,3):  u(2,1,3) >= -0
MCLin7(2,1,3):  -47 y(1,2,3) - u(2,1,3) <= -0
MCLin8(2,1,3):  36 y(1,2,3) + v(2,1,3) - u(2,1,3) <= 36
MCLin9(2,1,3):  36 y(1,2,3) - u(2,1,3) >= -0
MCLin10(2,1,3):  -47 y(1,2,3) + v(2,1,3) - u(2,1,3) >= -47
makespan(0):  C - x(0,0) >= 4
makespan(1):  C - x(1,0) >= 9
makespan(2):  C - x(2,0) >= 7
Bounds
 0 <= C <= 62
 17 <= x(0,0) <= 58
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -45 <= v(0,1,0) <= 30
 -43 <= v(1,0,0) <= 32
 -620 <= u(0,1,0)
 -620 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -45 <= v(0,2,0) <= 45
 -56 <= v(2,0,0) <= 34
 -620 <= u(0,2,0)
 -620 <= u(2,0,0)
 19 <= x(1,0) <= 53
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -43 <= v(1,2,0) <= 40
 -56 <= v(2,1,0) <= 27
 -620 <= u(1,2,0)
 -620 <= u(2,1,0)
 6 <= x(2,0) <= 55
 0 <= x(0,1) <= 41
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -42 <= v(0,1,1) <= 33
 -46 <= v(1,0,1) <= 29
 -620 <= u(0,1,1)
 -620 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -51 <= v(0,2,1) <= 39
 -46 <= v(2,0,1) <= 44
 -620 <= u(0,2,1)
 -620 <= u(2,0,1)
 0 <= x(1,1) <= 34
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -51 <= v(1,2,1) <= 32
 -42 <= v(2,1,1) <= 41
 -620 <= u(1,2,1)
 -620 <= u(2,1,1)
 0 <= x(2,1) <= 49
 5 <= x(0,2) <= 46
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -48 <= v(0,1,2) <= 27
 -30 <= v(1,0,2) <= 45
 -620 <= u(0,1,2)
 -620 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -49 <= v(0,2,2) <= 41
 -46 <= v(2,0,2) <= 44
 -620 <= u(0,2,2)
 -620 <= u(2,0,2)
 18 <= x(1,2) <= 52
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -36 <= v(1,2,2) <= 47
 -51 <= v(2,1,2) <= 32
 -620 <= u(1,2,2)
 -620 <= u(2,1,2)
 2 <= x(2,2) <= 51
 7 <= x(0,3) <= 48
 0 <= y(0,1,3) <= 1
 0 <= y(1,0,3) <= 1
 -45 <= v(0,1,3) <= 30
 -50 <= v(1,0,3) <= 25
 -620 <= u(0,1,3)
 -620 <= u(1,0,3)
 0 <= y(0,2,3) <= 1
 0 <= y(2,0,3) <= 1
 -48 <= v(0,2,3) <= 42
 -53 <= v(2,0,3) <= 37
 -620 <= u(0,2,3)
 -620 <= u(2,0,3)
 8 <= x(1,3) <= 42
 0 <= y(1,2,3) <= 1
 0 <= y(2,1,3) <= 1
 -47 <= v(1,2,3) <= 36
 -47 <= v(2,1,3) <= 36
 -620 <= u(1,2,3)
 -620 <= u(2,1,3)
 5 <= x(2,3) <= 54
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
