\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,0) + x(0,2) >= 3
ord(0,2):  x(0,1) - x(0,2) >= 4
ord(1,1):  x(1,0) - x(1,1) >= 7
ord(1,2):  - x(1,0) + x(1,2) >= 8
ord(2,1):  x(2,0) - x(2,1) >= 2
ord(2,2):  - x(2,0) + x(2,2) >= 6
ord(3,1):  x(3,0) - x(3,2) >= 7
ord(3,2):  - x(3,0) + x(3,1) >= 7
ord(4,1):  - x(4,1) + x(4,2) >= 5
ord(4,2):  x(4,0) - x(4,2) >= 4
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(0,3,0):  y(0,3,0) + y(3,0,0) = 1
triangle2(0,3,1):  y(0,3,1) + y(3,0,1) = 1
triangle2(0,3,2):  y(0,3,2) + y(3,0,2) = 1
triangle2(0,4,0):  y(0,4,0) + y(4,0,0) = 1
triangle2(0,4,1):  y(0,4,1) + y(4,0,1) = 1
triangle2(0,4,2):  y(0,4,2) + y(4,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
triangle2(1,3,0):  y(1,3,0) + y(3,1,0) = 1
triangle2(1,3,1):  y(1,3,1) + y(3,1,1) = 1
triangle2(1,3,2):  y(1,3,2) + y(3,1,2) = 1
triangle2(1,4,0):  y(1,4,0) + y(4,1,0) = 1
triangle2(1,4,1):  y(1,4,1) + y(4,1,1) = 1
triangle2(1,4,2):  y(1,4,2) + y(4,1,2) = 1
triangle2(2,3,0):  y(2,3,0) + y(3,2,0) = 1
triangle2(2,3,1):  y(2,3,1) + y(3,2,1) = 1
triangle2(2,3,2):  y(2,3,2) + y(3,2,2) = 1
triangle2(2,4,0):  y(2,4,0) + y(4,2,0) = 1
triangle2(2,4,1):  y(2,4,1) + y(4,2,1) = 1
triangle2(2,4,2):  y(2,4,2) + y(4,2,2) = 1
triangle2(3,4,0):  y(3,4,0) + y(4,3,0) = 1
triangle2(3,4,1):  y(3,4,1) + y(4,3,1) = 1
triangle2(3,4,2):  y(3,4,2) + y(4,3,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 8
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 3
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -66 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  37 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 37
MCLin4(0,1,0):  37 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -66 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -66
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -48 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  55 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 55
MCLin9(1,0,0):  55 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -48 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -48
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 7
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 8
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -51 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  52 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 52
MCLin4(0,1,1):  52 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -51 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -51
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -67 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  36 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 36
MCLin9(1,0,1):  36 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -67 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -67
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 1
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 4
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -64 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  39 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 39
MCLin4(0,1,2):  39 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -64 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -64
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -44 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  59 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 59
MCLin9(1,0,2):  59 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -44 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -44
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 6
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 3
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -64 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  44 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 44
MCLin4(0,2,0):  44 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -64 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -64
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -53 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  55 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 55
MCLin9(2,0,0):  55 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -53 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -53
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 8
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -51 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  57 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 57
MCLin4(0,2,1):  57 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -51 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -51
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -67 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  41 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 41
MCLin9(2,0,1):  41 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -67 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -67
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 3
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 4
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -64 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  44 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 44
MCLin4(0,2,2):  44 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -64 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -64
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -51 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  57 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 57
MCLin9(2,0,2):  57 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -51 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -51
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 7
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 3
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -66 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  38 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 38
MCLin4(0,3,0):  38 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -66 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -66
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -48 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  56 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 56
MCLin9(3,0,0):  56 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -48 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -48
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 1
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 8
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -60 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  44 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 44
MCLin4(0,3,1):  44 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -60 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -60
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -53 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  51 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 51
MCLin9(3,0,1):  51 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -53 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -53
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 7
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 4
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -56 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  48 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 48
MCLin4(0,3,2):  48 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -56 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -56
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -59 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  45 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 45
MCLin9(3,0,2):  45 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -59 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -59
MCLinV(0,4,0):  x(0,0) - v(0,4,0) - x(4,0) = 1
MCLinV(4,0,0):  - x(0,0) - v(4,0,0) + x(4,0) = 3
MCLin1(0,4,0):  u(0,4,0) >= -0
MCLin2(0,4,0):  -67 y(4,0,0) - u(0,4,0) <= -0
MCLin3(0,4,0):  42 y(4,0,0) + v(0,4,0) - u(0,4,0) <= 42
MCLin4(0,4,0):  42 y(4,0,0) - u(0,4,0) >= -0
MCLin5(0,4,0):  -67 y(4,0,0) + v(0,4,0) - u(0,4,0) >= -67
MCLin6(4,0,0):  u(4,0,0) >= -0
MCLin7(4,0,0):  -46 y(0,4,0) - u(4,0,0) <= -0
MCLin8(4,0,0):  63 y(0,4,0) + v(4,0,0) - u(4,0,0) <= 63
MCLin9(4,0,0):  63 y(0,4,0) - u(4,0,0) >= -0
MCLin10(4,0,0):  -46 y(0,4,0) + v(4,0,0) - u(4,0,0) >= -46
MCLinV(0,4,1):  x(0,1) - v(0,4,1) - x(4,1) = 5
MCLinV(4,0,1):  - x(0,1) - v(4,0,1) + x(4,1) = 8
MCLin1(0,4,1):  u(0,4,1) >= -0
MCLin2(0,4,1):  -55 y(4,0,1) - u(0,4,1) <= -0
MCLin3(0,4,1):  54 y(4,0,1) + v(0,4,1) - u(0,4,1) <= 54
MCLin4(0,4,1):  54 y(4,0,1) - u(0,4,1) >= -0
MCLin5(0,4,1):  -55 y(4,0,1) + v(0,4,1) - u(0,4,1) >= -55
MCLin6(4,0,1):  u(4,0,1) >= -0
MCLin7(4,0,1):  -67 y(0,4,1) - u(4,0,1) <= -0
MCLin8(4,0,1):  42 y(0,4,1) + v(4,0,1) - u(4,0,1) <= 42
MCLin9(4,0,1):  42 y(0,4,1) - u(4,0,1) >= -0
MCLin10(4,0,1):  -67 y(0,4,1) + v(4,0,1) - u(4,0,1) >= -67
MCLinV(0,4,2):  x(0,2) - v(0,4,2) - x(4,2) = 4
MCLinV(4,0,2):  - x(0,2) - v(4,0,2) + x(4,2) = 4
MCLin1(0,4,2):  u(0,4,2) >= -0
MCLin2(0,4,2):  -63 y(4,0,2) - u(0,4,2) <= -0
MCLin3(0,4,2):  46 y(4,0,2) + v(0,4,2) - u(0,4,2) <= 46
MCLin4(0,4,2):  46 y(4,0,2) - u(0,4,2) >= -0
MCLin5(0,4,2):  -63 y(4,0,2) + v(0,4,2) - u(0,4,2) >= -63
MCLin6(4,0,2):  u(4,0,2) >= -0
MCLin7(4,0,2):  -54 y(0,4,2) - u(4,0,2) <= -0
MCLin8(4,0,2):  55 y(0,4,2) + v(4,0,2) - u(4,0,2) <= 55
MCLin9(4,0,2):  55 y(0,4,2) - u(4,0,2) >= -0
MCLin10(4,0,2):  -54 y(0,4,2) + v(4,0,2) - u(4,0,2) >= -54
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 6
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 8
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -57 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  50 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 50
MCLin4(1,2,0):  50 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -57 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -57
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -64 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  43 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 43
MCLin9(2,1,0):  43 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -64 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -64
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 7
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -58 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  49 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 49
MCLin4(1,2,1):  49 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -58 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -58
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -58 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  49 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 49
MCLin9(2,1,1):  49 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -58 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -58
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 3
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 1
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -52 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  55 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 55
MCLin4(1,2,2):  55 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -52 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -52
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -59 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  48 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 48
MCLin9(2,1,2):  48 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -59 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -59
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 7
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 8
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -59 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  44 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 44
MCLin4(1,3,0):  44 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -59 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -59
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -59 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  44 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 44
MCLin9(3,1,0):  44 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -59 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -59
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 1
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 7
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -67 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  36 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 36
MCLin4(1,3,1):  36 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -67 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -67
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -44 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  59 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 59
MCLin9(3,1,1):  59 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -44 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -44
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 7
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 1
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -44 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  59 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 59
MCLin4(1,3,2):  59 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -44 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -44
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -67 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  36 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 36
MCLin9(3,1,2):  36 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -67 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -67
MCLinV(1,4,0):  x(1,0) - v(1,4,0) - x(4,0) = 1
MCLinV(4,1,0):  - x(1,0) - v(4,1,0) + x(4,0) = 8
MCLin1(1,4,0):  u(1,4,0) >= -0
MCLin2(1,4,0):  -60 y(4,1,0) - u(1,4,0) <= -0
MCLin3(1,4,0):  48 y(4,1,0) + v(1,4,0) - u(1,4,0) <= 48
MCLin4(1,4,0):  48 y(4,1,0) - u(1,4,0) >= -0
MCLin5(1,4,0):  -60 y(4,1,0) + v(1,4,0) - u(1,4,0) >= -60
MCLin6(4,1,0):  u(4,1,0) >= -0
MCLin7(4,1,0):  -57 y(1,4,0) - u(4,1,0) <= -0
MCLin8(4,1,0):  51 y(1,4,0) + v(4,1,0) - u(4,1,0) <= 51
MCLin9(4,1,0):  51 y(1,4,0) - u(4,1,0) >= -0
MCLin10(4,1,0):  -57 y(1,4,0) + v(4,1,0) - u(4,1,0) >= -57
MCLinV(1,4,1):  x(1,1) - v(1,4,1) - x(4,1) = 5
MCLinV(4,1,1):  - x(1,1) - v(4,1,1) + x(4,1) = 7
MCLin1(1,4,1):  u(1,4,1) >= -0
MCLin2(1,4,1):  -62 y(4,1,1) - u(1,4,1) <= -0
MCLin3(1,4,1):  46 y(4,1,1) + v(1,4,1) - u(1,4,1) <= 46
MCLin4(1,4,1):  46 y(4,1,1) - u(1,4,1) >= -0
MCLin5(1,4,1):  -62 y(4,1,1) + v(1,4,1) - u(1,4,1) >= -62
MCLin6(4,1,1):  u(4,1,1) >= -0
MCLin7(4,1,1):  -58 y(1,4,1) - u(4,1,1) <= -0
MCLin8(4,1,1):  50 y(1,4,1) + v(4,1,1) - u(4,1,1) <= 50
MCLin9(4,1,1):  50 y(1,4,1) - u(4,1,1) >= -0
MCLin10(4,1,1):  -58 y(1,4,1) + v(4,1,1) - u(4,1,1) >= -58
MCLinV(1,4,2):  x(1,2) - v(1,4,2) - x(4,2) = 4
MCLinV(4,1,2):  - x(1,2) - v(4,1,2) + x(4,2) = 1
MCLin1(1,4,2):  u(1,4,2) >= -0
MCLin2(1,4,2):  -51 y(4,1,2) - u(1,4,2) <= -0
MCLin3(1,4,2):  57 y(4,1,2) + v(1,4,2) - u(1,4,2) <= 57
MCLin4(1,4,2):  57 y(4,1,2) - u(1,4,2) >= -0
MCLin5(1,4,2):  -51 y(4,1,2) + v(1,4,2) - u(1,4,2) >= -51
MCLin6(4,1,2):  u(4,1,2) >= -0
MCLin7(4,1,2):  -62 y(1,4,2) - u(4,1,2) <= -0
MCLin8(4,1,2):  46 y(1,4,2) + v(4,1,2) - u(4,1,2) <= 46
MCLin9(4,1,2):  46 y(1,4,2) - u(4,1,2) >= -0
MCLin10(4,1,2):  -62 y(1,4,2) + v(4,1,2) - u(4,1,2) >= -62
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 7
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 6
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -64 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  44 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 44
MCLin4(2,3,0):  44 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -64 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -64
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -57 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  51 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 51
MCLin9(3,2,0):  51 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -57 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -57
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 1
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 2
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -67 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  41 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 41
MCLin4(2,3,1):  41 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -67 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -67
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -44 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  64 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 64
MCLin9(3,2,1):  64 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -44 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -44
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 7
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 3
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -51 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  57 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 57
MCLin4(2,3,2):  57 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -51 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -51
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -67 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  41 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 41
MCLin9(3,2,2):  41 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -67 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -67
MCLinV(2,4,0):  x(2,0) - v(2,4,0) - x(4,0) = 1
MCLinV(4,2,0):  - x(2,0) - v(4,2,0) + x(4,0) = 6
MCLin1(2,4,0):  u(2,4,0) >= -0
MCLin2(2,4,0):  -65 y(4,2,0) - u(2,4,0) <= -0
MCLin3(2,4,0):  48 y(4,2,0) + v(2,4,0) - u(2,4,0) <= 48
MCLin4(2,4,0):  48 y(4,2,0) - u(2,4,0) >= -0
MCLin5(2,4,0):  -65 y(4,2,0) + v(2,4,0) - u(2,4,0) >= -65
MCLin6(4,2,0):  u(4,2,0) >= -0
MCLin7(4,2,0):  -55 y(2,4,0) - u(4,2,0) <= -0
MCLin8(4,2,0):  58 y(2,4,0) + v(4,2,0) - u(4,2,0) <= 58
MCLin9(4,2,0):  58 y(2,4,0) - u(4,2,0) >= -0
MCLin10(4,2,0):  -55 y(2,4,0) + v(4,2,0) - u(4,2,0) >= -55
MCLinV(2,4,1):  x(2,1) - v(2,4,1) - x(4,1) = 5
MCLinV(4,2,1):  - x(2,1) - v(4,2,1) + x(4,1) = 2
MCLin1(2,4,1):  u(2,4,1) >= -0
MCLin2(2,4,1):  -62 y(4,2,1) - u(2,4,1) <= -0
MCLin3(2,4,1):  51 y(4,2,1) + v(2,4,1) - u(2,4,1) <= 51
MCLin4(2,4,1):  51 y(4,2,1) - u(2,4,1) >= -0
MCLin5(2,4,1):  -62 y(4,2,1) + v(2,4,1) - u(2,4,1) >= -62
MCLin6(4,2,1):  u(4,2,1) >= -0
MCLin7(4,2,1):  -58 y(2,4,1) - u(4,2,1) <= -0
MCLin8(4,2,1):  55 y(2,4,1) + v(4,2,1) - u(4,2,1) <= 55
MCLin9(4,2,1):  55 y(2,4,1) - u(4,2,1) >= -0
MCLin10(4,2,1):  -58 y(2,4,1) + v(4,2,1) - u(4,2,1) >= -58
MCLinV(2,4,2):  x(2,2) - v(2,4,2) - x(4,2) = 4
MCLinV(4,2,2):  - x(2,2) - v(4,2,2) + x(4,2) = 3
MCLin1(2,4,2):  u(2,4,2) >= -0
MCLin2(2,4,2):  -58 y(4,2,2) - u(2,4,2) <= -0
MCLin3(2,4,2):  55 y(4,2,2) + v(2,4,2) - u(2,4,2) <= 55
MCLin4(2,4,2):  55 y(4,2,2) - u(2,4,2) >= -0
MCLin5(2,4,2):  -58 y(4,2,2) + v(2,4,2) - u(2,4,2) >= -58
MCLin6(4,2,2):  u(4,2,2) >= -0
MCLin7(4,2,2):  -62 y(2,4,2) - u(4,2,2) <= -0
MCLin8(4,2,2):  51 y(2,4,2) + v(4,2,2) - u(4,2,2) <= 51
MCLin9(4,2,2):  51 y(2,4,2) - u(4,2,2) >= -0
MCLin10(4,2,2):  -62 y(2,4,2) + v(4,2,2) - u(4,2,2) >= -62
MCLinV(3,4,0):  x(3,0) - v(3,4,0) - x(4,0) = 1
MCLinV(4,3,0):  - x(3,0) - v(4,3,0) + x(4,0) = 7
MCLin1(3,4,0):  u(3,4,0) >= -0
MCLin2(3,4,0):  -60 y(4,3,0) - u(3,4,0) <= -0
MCLin3(3,4,0):  49 y(4,3,0) + v(3,4,0) - u(3,4,0) <= 49
MCLin4(3,4,0):  49 y(4,3,0) - u(3,4,0) >= -0
MCLin5(3,4,0):  -60 y(4,3,0) + v(3,4,0) - u(3,4,0) >= -60
MCLin6(4,3,0):  u(4,3,0) >= -0
MCLin7(4,3,0):  -57 y(3,4,0) - u(4,3,0) <= -0
MCLin8(4,3,0):  52 y(3,4,0) + v(4,3,0) - u(4,3,0) <= 52
MCLin9(4,3,0):  52 y(3,4,0) - u(4,3,0) >= -0
MCLin10(4,3,0):  -57 y(3,4,0) + v(4,3,0) - u(4,3,0) >= -57
MCLinV(3,4,1):  x(3,1) - v(3,4,1) - x(4,1) = 5
MCLinV(4,3,1):  - x(3,1) - v(4,3,1) + x(4,1) = 1
MCLin1(3,4,1):  u(3,4,1) >= -0
MCLin2(3,4,1):  -48 y(4,3,1) - u(3,4,1) <= -0
MCLin3(3,4,1):  61 y(4,3,1) + v(3,4,1) - u(3,4,1) <= 61
MCLin4(3,4,1):  61 y(4,3,1) - u(3,4,1) >= -0
MCLin5(3,4,1):  -48 y(4,3,1) + v(3,4,1) - u(3,4,1) >= -48
MCLin6(4,3,1):  u(4,3,1) >= -0
MCLin7(4,3,1):  -67 y(3,4,1) - u(4,3,1) <= -0
MCLin8(4,3,1):  42 y(3,4,1) + v(4,3,1) - u(4,3,1) <= 42
MCLin9(4,3,1):  42 y(3,4,1) - u(4,3,1) >= -0
MCLin10(4,3,1):  -67 y(3,4,1) + v(4,3,1) - u(4,3,1) >= -67
MCLinV(3,4,2):  x(3,2) - v(3,4,2) - x(4,2) = 4
MCLinV(4,3,2):  - x(3,2) - v(4,3,2) + x(4,2) = 7
MCLin1(3,4,2):  u(3,4,2) >= -0
MCLin2(3,4,2):  -66 y(4,3,2) - u(3,4,2) <= -0
MCLin3(3,4,2):  43 y(4,3,2) + v(3,4,2) - u(3,4,2) <= 43
MCLin4(3,4,2):  43 y(4,3,2) - u(3,4,2) >= -0
MCLin5(3,4,2):  -66 y(4,3,2) + v(3,4,2) - u(3,4,2) >= -66
MCLin6(4,3,2):  u(4,3,2) >= -0
MCLin7(4,3,2):  -54 y(3,4,2) - u(4,3,2) <= -0
MCLin8(4,3,2):  55 y(3,4,2) + v(4,3,2) - u(4,3,2) <= 55
MCLin9(4,3,2):  55 y(3,4,2) - u(4,3,2) >= -0
MCLin10(4,3,2):  -54 y(3,4,2) + v(4,3,2) - u(4,3,2) >= -54
makespan(0):  C - x(0,1) >= 8
makespan(1):  C - x(1,2) >= 1
makespan(2):  C - x(2,2) >= 3
makespan(3):  C - x(3,1) >= 1
makespan(4):  C - x(4,0) >= 1
cut_clique_0(01234,0):  0.00334 x(0,0) + 0.03640 x(1,0) + 0.00728 x(2,0) + 0.03185 x(3,0) + 0.00455 x(4,0) >= 1
cut_clique_0(1234,0):  0.03441 x(1,0) + 0.00369 x(2,0) + 0.02904 x(3,0) + 0.01398 x(4,0) >= 1
cut_clique_0(0123,0):  0.00368 x(0,0) + 0.04087 x(1,0) + 0.00807 x(2,0) + 0.03576 x(3,0) >= 1
cut_clique_0(0134,0):  0.03187 x(1,0) + 0.02561 x(3,0) + 0.02447 x(4,0) >= 1
cut_clique_0(0124,0):  0.00656 x(0,0) + 0.08588 x(1,0) + 0.01486 x(2,0) + 0.01073 x(4,0) >= 1
cut_clique_0(0234,0):  0.00703 x(0,0) + 0.01606 x(2,0) + 0.08434 x(3,0) + 0.01205 x(4,0) >= 1
cut_clique_0(134,0):  0.03187 x(1,0) + 0.02561 x(3,0) + 0.02447 x(4,0) >= 1
cut_clique_0(123,0):  0.04507 x(1,0) + 0.00423 x(2,0) + 0.03944 x(3,0) >= 1
cut_clique_0(124,0):  0.07840 x(1,0) + 0.00720 x(2,0) + 0.02240 x(4,0) >= 1
cut_clique_0(234,0):  0.00777 x(2,0) + 0.07573 x(3,0) + 0.02524 x(4,0) >= 1
cut_clique_0(012,0):  0.00768 x(0,0) + 0.10526 x(1,0) + 0.01754 x(2,0) >= 1
cut_clique_0(013,0):  0.04969 x(1,0) + 0.04348 x(3,0) >= 1
cut_clique_0(023,0):  0.00825 x(0,0) + 0.01905 x(2,0) + 0.10476 x(3,0) >= 1
cut_clique_0(014,0):  0.06897 x(1,0) + 0.03448 x(4,0) >= 1
cut_clique_0(034,0):  0.06494 x(3,0) + 0.03896 x(4,0) >= 1
cut_clique_0(01234,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(0124,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(1234,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(0123,1):  0.50000 x(1,1) + 0.14286 x(2,1) >= 1
cut_clique_0(0234,1):  0.20000 x(2,1) + 0.50000 x(4,1) >= 1
cut_clique_0(0134,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(124,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(012,1):  0.50000 x(1,1) + 0.14286 x(2,1) >= 1
cut_clique_0(014,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(024,1):  0.20000 x(2,1) + 0.50000 x(4,1) >= 1
cut_clique_0(123,1):  0.50000 x(1,1) + 0.14286 x(2,1) >= 1
cut_clique_0(134,1):  0.20000 x(1,1) + 0.14286 x(4,1) >= 1
cut_clique_0(234,1):  0.20000 x(2,1) + 0.50000 x(4,1) >= 1
cut_clique_0(01234,2):  0.10116 x(0,2) + 0.04643 x(3,2) + 0.02653 x(4,2) >= 1
cut_clique_0(0234,2):  0.10116 x(0,2) + 0.04643 x(3,2) + 0.02653 x(4,2) >= 1
cut_clique_0(0124,2):  0.04167 x(0,2) + 0.12500 x(4,2) >= 1
cut_clique_0(0134,2):  0.10116 x(0,2) + 0.04643 x(3,2) + 0.02653 x(4,2) >= 1
cut_clique_0(1234,2):  0.03175 x(3,2) + 0.14286 x(4,2) >= 1
cut_clique_0(0123,2):  0.14286 x(0,2) + 0.08163 x(3,2) >= 1
cut_clique_0(034,2):  0.10116 x(0,2) + 0.04643 x(3,2) + 0.02653 x(4,2) >= 1
cut_clique_0(024,2):  0.04167 x(0,2) + 0.12500 x(4,2) >= 1
cut_clique_0(023,2):  0.14286 x(0,2) + 0.08163 x(3,2) >= 1
cut_clique_0(234,2):  0.03175 x(3,2) + 0.14286 x(4,2) >= 1
cut_clique_0(014,2):  0.04167 x(0,2) + 0.12500 x(4,2) >= 1
cut_clique_0(013,2):  0.14286 x(0,2) + 0.08163 x(3,2) >= 1
cut_clique_0(134,2):  0.03175 x(3,2) + 0.14286 x(4,2) >= 1
cut_clique_1(01234,1):  0.04145 x(0,1) + 0.03627 x(1,1) + 0.01036 x(2,1) + 0.00518 x(3,1) + 0.02591 x(4,1) >= 1
cut_clique_1(0124,1):  0.04678 x(0,1) + 0.04094 x(1,1) + 0.01170 x(2,1) + 0.02924 x(4,1) >= 1
cut_clique_1(1234,1):  0.11864 x(1,1) + 0.03390 x(2,1) + 0.08475 x(4,1) >= 1
cut_clique_1(124,1):  0.11864 x(1,1) + 0.03390 x(2,1) + 0.08475 x(4,1) >= 1
cut_clique_1(01234,2):  0.05556 x(0,2) + 0.04040 x(3,2) + 0.05556 x(4,2) >= 1
cut_clique_1(0234,2):  0.05556 x(0,2) + 0.04040 x(3,2) + 0.05556 x(4,2) >= 1
cut_clique_1(0134,2):  0.05556 x(0,2) + 0.04040 x(3,2) + 0.05556 x(4,2) >= 1
cut_clique_1(034,2):  0.05556 x(0,2) + 0.04040 x(3,2) + 0.05556 x(4,2) >= 1
cut_clique_0:  0.03276 x(1,0) + 0.00655 x(2,0) + 0.02867 x(3,0) + 0.00410 x(4,0) + 0.00554 x(0,2) + 0.00403 x(3,2) + 0.00554 x(4,2) >= 1
cut_clique_1:  0.03523 x(1,0) + 0.00705 x(2,0) + 0.03083 x(3,0) + 0.00440 x(4,0) + 0.00456 x(0,2) + 0.00261 x(3,2) >= 1
Bounds
 0 <= C <= 67
 0 <= x(0,0) <= 52
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -66 <= v(0,1,0) <= 37
 -48 <= v(1,0,0) <= 55
 -670 <= u(0,1,0)
 -670 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -64 <= v(0,2,0) <= 44
 -53 <= v(2,0,0) <= 55
 -670 <= u(0,2,0)
 -670 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -66 <= v(0,3,0) <= 38
 -48 <= v(3,0,0) <= 56
 -670 <= u(0,3,0)
 -670 <= u(3,0,0)
 0 <= y(0,4,0) <= 1
 0 <= y(4,0,0) <= 1
 -67 <= v(0,4,0) <= 42
 -46 <= v(4,0,0) <= 63
 -670 <= u(0,4,0)
 -670 <= u(4,0,0)
 7 <= x(1,0) <= 58
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -57 <= v(1,2,0) <= 50
 -64 <= v(2,1,0) <= 43
 -670 <= u(1,2,0)
 -670 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -59 <= v(1,3,0) <= 44
 -59 <= v(3,1,0) <= 44
 -670 <= u(1,3,0)
 -670 <= u(3,1,0)
 0 <= y(1,4,0) <= 1
 0 <= y(4,1,0) <= 1
 -60 <= v(1,4,0) <= 48
 -57 <= v(4,1,0) <= 51
 -670 <= u(1,4,0)
 -670 <= u(4,1,0)
 2 <= x(2,0) <= 58
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -64 <= v(2,3,0) <= 44
 -57 <= v(3,2,0) <= 51
 -670 <= u(2,3,0)
 -670 <= u(3,2,0)
 0 <= y(2,4,0) <= 1
 0 <= y(4,2,0) <= 1
 -65 <= v(2,4,0) <= 48
 -55 <= v(4,2,0) <= 58
 -670 <= u(2,4,0)
 -670 <= u(4,2,0)
 7 <= x(3,0) <= 59
 0 <= y(3,4,0) <= 1
 0 <= y(4,3,0) <= 1
 -60 <= v(3,4,0) <= 49
 -57 <= v(4,3,0) <= 52
 -670 <= u(3,4,0)
 -670 <= u(4,3,0)
 9 <= x(4,0) <= 66
 7 <= x(0,1) <= 59
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -51 <= v(0,1,1) <= 52
 -67 <= v(1,0,1) <= 36
 -670 <= u(0,1,1)
 -670 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -51 <= v(0,2,1) <= 57
 -67 <= v(2,0,1) <= 41
 -670 <= u(0,2,1)
 -670 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -60 <= v(0,3,1) <= 44
 -53 <= v(3,0,1) <= 51
 -670 <= u(0,3,1)
 -670 <= u(3,0,1)
 0 <= y(0,4,1) <= 1
 0 <= y(4,0,1) <= 1
 -55 <= v(0,4,1) <= 54
 -67 <= v(4,0,1) <= 42
 -670 <= u(0,4,1)
 -670 <= u(4,0,1)
 0 <= x(1,1) <= 51
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -58 <= v(1,2,1) <= 49
 -58 <= v(2,1,1) <= 49
 -670 <= u(1,2,1)
 -670 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -67 <= v(1,3,1) <= 36
 -44 <= v(3,1,1) <= 59
 -670 <= u(1,3,1)
 -670 <= u(3,1,1)
 0 <= y(1,4,1) <= 1
 0 <= y(4,1,1) <= 1
 -62 <= v(1,4,1) <= 46
 -58 <= v(4,1,1) <= 50
 -670 <= u(1,4,1)
 -670 <= u(4,1,1)
 0 <= x(2,1) <= 56
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -67 <= v(2,3,1) <= 41
 -44 <= v(3,2,1) <= 64
 -670 <= u(2,3,1)
 -670 <= u(3,2,1)
 0 <= y(2,4,1) <= 1
 0 <= y(4,2,1) <= 1
 -62 <= v(2,4,1) <= 51
 -58 <= v(4,2,1) <= 55
 -670 <= u(2,4,1)
 -670 <= u(4,2,1)
 14 <= x(3,1) <= 66
 0 <= y(3,4,1) <= 1
 0 <= y(4,3,1) <= 1
 -48 <= v(3,4,1) <= 61
 -67 <= v(4,3,1) <= 42
 -670 <= u(3,4,1)
 -670 <= u(4,3,1)
 0 <= x(4,1) <= 57
 3 <= x(0,2) <= 55
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -64 <= v(0,1,2) <= 39
 -44 <= v(1,0,2) <= 59
 -670 <= u(0,1,2)
 -670 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -64 <= v(0,2,2) <= 44
 -51 <= v(2,0,2) <= 57
 -670 <= u(0,2,2)
 -670 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -56 <= v(0,3,2) <= 48
 -59 <= v(3,0,2) <= 45
 -670 <= u(0,3,2)
 -670 <= u(3,0,2)
 0 <= y(0,4,2) <= 1
 0 <= y(4,0,2) <= 1
 -63 <= v(0,4,2) <= 46
 -54 <= v(4,0,2) <= 55
 -670 <= u(0,4,2)
 -670 <= u(4,0,2)
 15 <= x(1,2) <= 66
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -52 <= v(1,2,2) <= 55
 -59 <= v(2,1,2) <= 48
 -670 <= u(1,2,2)
 -670 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -44 <= v(1,3,2) <= 59
 -67 <= v(3,1,2) <= 36
 -670 <= u(1,3,2)
 -670 <= u(3,1,2)
 0 <= y(1,4,2) <= 1
 0 <= y(4,1,2) <= 1
 -51 <= v(1,4,2) <= 57
 -62 <= v(4,1,2) <= 46
 -670 <= u(1,4,2)
 -670 <= u(4,1,2)
 8 <= x(2,2) <= 64
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -51 <= v(2,3,2) <= 57
 -67 <= v(3,2,2) <= 41
 -670 <= u(2,3,2)
 -670 <= u(3,2,2)
 0 <= y(2,4,2) <= 1
 0 <= y(4,2,2) <= 1
 -58 <= v(2,4,2) <= 55
 -62 <= v(4,2,2) <= 51
 -670 <= u(2,4,2)
 -670 <= u(4,2,2)
 0 <= x(3,2) <= 52
 0 <= y(3,4,2) <= 1
 0 <= y(4,3,2) <= 1
 -66 <= v(3,4,2) <= 43
 -54 <= v(4,3,2) <= 55
 -670 <= u(3,4,2)
 -670 <= u(4,3,2)
 5 <= x(4,2) <= 62
Integers
C x(0,0) y(0,1,0) y(1,0,0) v(0,1,0) v(1,0,0) u(0,1,0) u(1,0,0) y(0,2,0) y(2,0,0) 
v(0,2,0) v(2,0,0) u(0,2,0) u(2,0,0) y(0,3,0) y(3,0,0) v(0,3,0) v(3,0,0) u(0,3,0) u(3,0,0) 
y(0,4,0) y(4,0,0) v(0,4,0) v(4,0,0) u(0,4,0) u(4,0,0) x(1,0) y(1,2,0) y(2,1,0) v(1,2,0) 
v(2,1,0) u(1,2,0) u(2,1,0) y(1,3,0) y(3,1,0) v(1,3,0) v(3,1,0) u(1,3,0) u(3,1,0) y(1,4,0) 
y(4,1,0) v(1,4,0) v(4,1,0) u(1,4,0) u(4,1,0) x(2,0) y(2,3,0) y(3,2,0) v(2,3,0) v(3,2,0) 
u(2,3,0) u(3,2,0) y(2,4,0) y(4,2,0) v(2,4,0) v(4,2,0) u(2,4,0) u(4,2,0) x(3,0) y(3,4,0) 
y(4,3,0) v(3,4,0) v(4,3,0) u(3,4,0) u(4,3,0) x(4,0) x(0,1) y(0,1,1) y(1,0,1) v(0,1,1) 
v(1,0,1) u(0,1,1) u(1,0,1) y(0,2,1) y(2,0,1) v(0,2,1) v(2,0,1) u(0,2,1) u(2,0,1) y(0,3,1) 
y(3,0,1) v(0,3,1) v(3,0,1) u(0,3,1) u(3,0,1) y(0,4,1) y(4,0,1) v(0,4,1) v(4,0,1) u(0,4,1) 
u(4,0,1) x(1,1) y(1,2,1) y(2,1,1) v(1,2,1) v(2,1,1) u(1,2,1) u(2,1,1) y(1,3,1) y(3,1,1) 
v(1,3,1) v(3,1,1) u(1,3,1) u(3,1,1) y(1,4,1) y(4,1,1) v(1,4,1) v(4,1,1) u(1,4,1) u(4,1,1) 
x(2,1) y(2,3,1) y(3,2,1) v(2,3,1) v(3,2,1) u(2,3,1) u(3,2,1) y(2,4,1) y(4,2,1) v(2,4,1) 
v(4,2,1) u(2,4,1) u(4,2,1) x(3,1) y(3,4,1) y(4,3,1) v(3,4,1) v(4,3,1) u(3,4,1) u(4,3,1) 
x(4,1) x(0,2) y(0,1,2) y(1,0,2) v(0,1,2) v(1,0,2) u(0,1,2) u(1,0,2) y(0,2,2) y(2,0,2) 
v(0,2,2) v(2,0,2) u(0,2,2) u(2,0,2) y(0,3,2) y(3,0,2) v(0,3,2) v(3,0,2) u(0,3,2) u(3,0,2) 
y(0,4,2) y(4,0,2) v(0,4,2) v(4,0,2) u(0,4,2) u(4,0,2) x(1,2) y(1,2,2) y(2,1,2) v(1,2,2) 
v(2,1,2) u(1,2,2) u(2,1,2) y(1,3,2) y(3,1,2) v(1,3,2) v(3,1,2) u(1,3,2) u(3,1,2) y(1,4,2) 
y(4,1,2) v(1,4,2) v(4,1,2) u(1,4,2) u(4,1,2) x(2,2) y(2,3,2) y(3,2,2) v(2,3,2) v(3,2,2) 
u(2,3,2) u(3,2,2) y(2,4,2) y(4,2,2) v(2,4,2) v(4,2,2) u(2,4,2) u(4,2,2) x(3,2) y(3,4,2) 
y(4,3,2) v(3,4,2) v(4,3,2) u(3,4,2) u(4,3,2) x(4,2) 
End
