\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,2) >= 4
ord(0,2):  x(0,0) - x(0,2) >= 5
ord(1,1):  - x(1,1) + x(1,2) >= 8
ord(1,2):  x(1,0) - x(1,2) >= 4
ord(2,1):  - x(2,1) + x(2,2) >= 2
ord(2,2):  x(2,0) - x(2,2) >= 1
ord(3,1):  - x(3,0) + x(3,1) >= 6
ord(3,2):  - x(3,1) + x(3,2) >= 2
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(0,3,0):  y(0,3,0) + y(3,0,0) = 1
triangle2(0,3,1):  y(0,3,1) + y(3,0,1) = 1
triangle2(0,3,2):  y(0,3,2) + y(3,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
triangle2(1,3,0):  y(1,3,0) + y(3,1,0) = 1
triangle2(1,3,1):  y(1,3,1) + y(3,1,1) = 1
triangle2(1,3,2):  y(1,3,2) + y(3,1,2) = 1
triangle2(2,3,0):  y(2,3,0) + y(3,2,0) = 1
triangle2(2,3,1):  y(2,3,1) + y(3,2,1) = 1
triangle2(2,3,2):  y(2,3,2) + y(3,2,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 10
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 7
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -52 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  32 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 32
MCLin4(0,1,0):  32 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -52 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -52
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -49 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  35 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 35
MCLin9(1,0,0):  35 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -49 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -49
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 8
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 4
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -47 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  37 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 37
MCLin4(0,1,1):  37 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -47 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -47
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -49 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  35 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 35
MCLin9(1,0,1):  35 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -49 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -49
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 4
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 5
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -47 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  37 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 37
MCLin4(0,1,2):  37 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -47 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -47
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -46 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  38 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 38
MCLin9(1,0,2):  38 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -46 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -46
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 3
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 7
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -52 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  48 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 48
MCLin4(0,2,0):  48 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -52 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -52
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -58 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  42 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 42
MCLin9(2,0,0):  42 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -58 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -58
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 4
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -57 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  43 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 43
MCLin4(0,2,1):  43 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -57 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -57
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -49 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  51 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 51
MCLin9(2,0,1):  51 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -49 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -49
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 1
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 5
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -54 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  46 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 46
MCLin4(0,2,2):  46 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -54 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -54
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -52 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  48 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 48
MCLin9(2,0,2):  48 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -52 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -52
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 6
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 7
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -41 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  48 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 48
MCLin4(0,3,0):  48 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -41 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -41
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -61 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  28 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 28
MCLin9(3,0,0):  28 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -61 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -61
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 2
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 4
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -52 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  37 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 37
MCLin4(0,3,1):  37 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -52 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -52
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -43 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  46 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 46
MCLin9(3,0,1):  46 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -43 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -43
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 9
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 5
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -57 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  32 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 32
MCLin4(0,3,2):  32 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -57 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -57
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -46 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  43 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 43
MCLin9(3,0,2):  43 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -46 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -46
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 3
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 10
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -49 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  45 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 45
MCLin4(1,2,0):  45 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -49 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -49
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -58 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  36 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 36
MCLin9(2,1,0):  36 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -58 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -58
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 8
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -57 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  37 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 37
MCLin4(1,2,1):  37 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -57 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -57
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -47 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  47 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 47
MCLin9(2,1,1):  47 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -47 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -47
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 1
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 4
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -50 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  44 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 44
MCLin4(1,2,2):  44 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -50 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -50
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -49 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  45 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 45
MCLin9(2,1,2):  45 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -49 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -49
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 6
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 10
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -38 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  45 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 45
MCLin4(1,3,0):  45 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -38 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -38
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -61 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  22 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 22
MCLin9(3,1,0):  22 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -61 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -61
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 2
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 8
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -52 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  31 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 31
MCLin4(1,3,1):  31 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -52 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -52
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -41 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  42 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 42
MCLin9(3,1,1):  42 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -41 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -41
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 9
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 4
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -53 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  30 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 30
MCLin4(1,3,2):  30 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -53 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -53
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -43 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  40 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 40
MCLin9(3,1,2):  40 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -43 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -43
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 6
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 3
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -47 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  52 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 52
MCLin4(2,3,0):  52 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -47 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -47
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -61 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  38 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 38
MCLin9(3,2,0):  38 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -61 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -61
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 2
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 2
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -52 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  47 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 47
MCLin4(2,3,1):  47 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -52 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -52
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -51 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  48 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 48
MCLin9(3,2,1):  48 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -51 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -51
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 9
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 1
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -59 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  40 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 40
MCLin4(2,3,2):  40 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -59 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -59
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -50 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  49 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 49
MCLin9(3,2,2):  49 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -50 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -50
makespan(0):  C - x(0,0) >= 7
makespan(1):  C - x(1,0) >= 10
makespan(2):  C - x(2,0) >= 3
makespan(3):  C - x(3,2) >= 9
cut_clique_0(0123,0):  0.01639 x(0,0) + 0.05328 x(1,0) >= 1
cut_clique_0(012,0):  0.01639 x(0,0) + 0.05328 x(1,0) >= 1
cut_clique_0(013,0):  0.01639 x(0,0) + 0.05328 x(1,0) >= 1
cut_clique_0(0123,1):  0.12500 x(0,1) + 0.25000 x(1,1) >= 1
cut_clique_0(012,1):  0.12500 x(0,1) + 0.25000 x(1,1) >= 1
cut_clique_0(013,1):  0.12500 x(0,1) + 0.25000 x(1,1) >= 1
cut_clique_0(023,1):  0.50000 x(0,1) + 0.25000 x(2,1) >= 1
cut_clique_0(123,1):  0.50000 x(1,1) + 0.12500 x(2,1) >= 1
cut_clique_0(0123,2):  0.00490 x(0,2) + 0.02563 x(1,2) + 0.05767 x(3,2) >= 1
cut_clique_0(123,2):  0.02857 x(1,2) + 0.06429 x(3,2) >= 1
cut_clique_0(012,2):  0.01316 x(0,2) + 0.10526 x(1,2) >= 1
cut_clique_0(013,2):  0.00490 x(0,2) + 0.02563 x(1,2) + 0.05767 x(3,2) >= 1
cut_clique_0(023,2):  0.00826 x(0,2) + 0.10744 x(3,2) >= 1
cut_clique_1(0123,1):  0.04762 x(0,1) + 0.09524 x(1,1) + 0.02381 x(2,1) + 0.02381 x(3,1) >= 1
cut_clique_1(013,1):  0.05556 x(0,1) + 0.08333 x(1,1) + 0.05556 x(3,1) >= 1
cut_clique_1(012,1):  0.07143 x(0,1) + 0.14286 x(1,1) + 0.03571 x(2,1) >= 1
cut_clique_1(123,1):  0.03125 x(1,1) + 0.12500 x(3,1) >= 1
cut_clique_0:  0.01212 x(0,0) + 0.03939 x(1,0) + 0.00116 x(0,2) + 0.01970 x(3,2) >= 1
cut_clique_1:  0.07347 x(1,1) + 0.02041 x(0,2) + 0.01633 x(1,2) + 0.03673 x(3,2) >= 1
cut_clique_2:  0.07143 x(0,1) + 0.14286 x(1,1) + 0.03571 x(3,1) >= 1
Bounds
 0 <= C <= 61
 9 <= x(0,0) <= 54
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -52 <= v(0,1,0) <= 32
 -49 <= v(1,0,0) <= 35
 -610 <= u(0,1,0)
 -610 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -52 <= v(0,2,0) <= 48
 -58 <= v(2,0,0) <= 42
 -610 <= u(0,2,0)
 -610 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -41 <= v(0,3,0) <= 48
 -61 <= v(3,0,0) <= 28
 -610 <= u(0,3,0)
 -610 <= u(3,0,0)
 12 <= x(1,0) <= 51
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -49 <= v(1,2,0) <= 45
 -58 <= v(2,1,0) <= 36
 -610 <= u(1,2,0)
 -610 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -38 <= v(1,3,0) <= 45
 -61 <= v(3,1,0) <= 22
 -610 <= u(1,3,0)
 -610 <= u(3,1,0)
 3 <= x(2,0) <= 58
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -47 <= v(2,3,0) <= 52
 -61 <= v(3,2,0) <= 38
 -610 <= u(2,3,0)
 -610 <= u(3,2,0)
 0 <= x(3,0) <= 44
 0 <= x(0,1) <= 45
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -47 <= v(0,1,1) <= 37
 -49 <= v(1,0,1) <= 35
 -610 <= u(0,1,1)
 -610 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -57 <= v(0,2,1) <= 43
 -49 <= v(2,0,1) <= 51
 -610 <= u(0,2,1)
 -610 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -52 <= v(0,3,1) <= 37
 -43 <= v(3,0,1) <= 46
 -610 <= u(0,3,1)
 -610 <= u(3,0,1)
 0 <= x(1,1) <= 39
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -57 <= v(1,2,1) <= 37
 -47 <= v(2,1,1) <= 47
 -610 <= u(1,2,1)
 -610 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -52 <= v(1,3,1) <= 31
 -41 <= v(3,1,1) <= 42
 -610 <= u(1,3,1)
 -610 <= u(3,1,1)
 0 <= x(2,1) <= 55
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -52 <= v(2,3,1) <= 47
 -51 <= v(3,2,1) <= 48
 -610 <= u(2,3,1)
 -610 <= u(3,2,1)
 6 <= x(3,1) <= 50
 4 <= x(0,2) <= 49
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -47 <= v(0,1,2) <= 37
 -46 <= v(1,0,2) <= 38
 -610 <= u(0,1,2)
 -610 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -54 <= v(0,2,2) <= 46
 -52 <= v(2,0,2) <= 48
 -610 <= u(0,2,2)
 -610 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -57 <= v(0,3,2) <= 32
 -46 <= v(3,0,2) <= 43
 -610 <= u(0,3,2)
 -610 <= u(3,0,2)
 8 <= x(1,2) <= 47
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -50 <= v(1,2,2) <= 44
 -49 <= v(2,1,2) <= 45
 -610 <= u(1,2,2)
 -610 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -53 <= v(1,3,2) <= 30
 -43 <= v(3,1,2) <= 40
 -610 <= u(1,3,2)
 -610 <= u(3,1,2)
 2 <= x(2,2) <= 57
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -59 <= v(2,3,2) <= 40
 -50 <= v(3,2,2) <= 49
 -610 <= u(2,3,2)
 -610 <= u(3,2,2)
 8 <= x(3,2) <= 52
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) y(0,3,0) y(3,0,0) v(0,3,0) v(3,0,0) u(0,3,0) u(3,0,0) 
x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) v(2,1,0) u(1,2,0) u(2,1,0) y(1,3,0) y(3,1,0) v(1,3,0) 
v(3,1,0) u(1,3,0) u(3,1,0) x(2,0) y(2,3,0) y(3,2,0) v(2,3,0) v(3,2,0) u(2,3,0) u(3,2,0) 
x(3,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) y(2,0,1) 
v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) y(0,3,1) y(3,0,1) v(0,3,1) v(3,0,1) u(0,3,1) u(3,0,1) 
x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) u(1,2,1) u(2,1,1) y(1,3,1) y(3,1,1) v(1,3,1) 
v(3,1,1) u(1,3,1) u(3,1,1) x(2,1) y(2,3,1) y(3,2,1) v(2,3,1) v(3,2,1) u(2,3,1) u(3,2,1) 
x(3,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) y(0,2,2) y(2,0,2) 
v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) y(0,3,2) y(3,0,2) v(0,3,2) v(3,0,2) u(0,3,2) u(3,0,2) 
x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) v(2,1,2) u(1,2,2) u(2,1,2) y(1,3,2) y(3,1,2) v(1,3,2) 
v(3,1,2) u(1,3,2) u(3,1,2) x(2,2) y(2,3,2) y(3,2,2) v(2,3,2) v(3,2,2) u(2,3,2) u(3,2,2) 
x(3,2) 
End
