\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  x(0,1) - x(0,2) >= 6
ord(0,2):  x(0,0) - x(0,1) >= 5
ord(1,1):  - x(1,0) + x(1,2) >= 4
ord(1,2):  x(1,1) - x(1,2) >= 4
ord(2,1):  - x(2,1) + x(2,2) >= 8
ord(2,2):  x(2,0) - x(2,2) >= 3
ord(3,1):  - x(3,1) + x(3,2) >= 4
ord(3,2):  x(3,0) - x(3,2) >= 10
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
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 4
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 1
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -32 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  43 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 43
MCLin4(0,1,0):  43 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -32 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -32
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -48 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  27 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 27
MCLin9(1,0,0):  27 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -48 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -48
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 1
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 5
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -42 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  33 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 33
MCLin4(0,1,1):  33 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -42 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -42
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -39 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  36 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 36
MCLin9(1,0,1):  36 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -39 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -39
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 4
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 6
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -47 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  28 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 28
MCLin4(0,1,2):  28 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -47 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -47
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -38 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  37 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 37
MCLin9(1,0,2):  37 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -38 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -38
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 1
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 1
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -37 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  35 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 35
MCLin4(0,2,0):  35 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -37 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -37
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -37 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  35 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 35
MCLin9(2,0,0):  35 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -37 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -37
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 8
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 5
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -38 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  34 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 34
MCLin4(0,2,1):  34 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -38 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -38
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -47 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  25 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 25
MCLin9(2,0,1):  25 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -47 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -47
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 3
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 6
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -47 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  25 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 25
MCLin4(0,2,2):  25 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -47 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -47
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -34 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  38 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 38
MCLin9(2,0,2):  38 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -34 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -34
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 1
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 1
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -37 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  32 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 32
MCLin4(0,3,0):  32 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -37 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -37
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -34 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  35 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 35
MCLin9(3,0,0):  35 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -34 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -34
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 4
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 5
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -31 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  38 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 38
MCLin4(0,3,1):  38 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -31 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -31
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -47 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  22 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 22
MCLin9(3,0,1):  22 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -47 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -47
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 10
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 6
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -47 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  22 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 22
MCLin4(0,3,2):  22 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -47 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -47
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -38 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  31 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 31
MCLin9(3,0,2):  31 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -38 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -38
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 1
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 4
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -48 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  27 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 27
MCLin4(1,2,0):  27 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -48 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -48
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -32 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  43 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 43
MCLin9(2,1,0):  43 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -32 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -32
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 8
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 1
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -36 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  39 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 39
MCLin4(1,2,1):  39 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -36 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -36
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -48 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  27 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 27
MCLin9(2,1,1):  27 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -48 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -48
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 3
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 4
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -43 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  32 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 32
MCLin4(1,2,2):  32 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -43 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -43
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -39 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  36 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 36
MCLin9(2,1,2):  36 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -39 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -39
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 1
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 4
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -48 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  24 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 24
MCLin4(1,3,0):  24 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -48 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -48
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -29 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  43 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 43
MCLin9(3,1,0):  43 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -29 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -29
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 4
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 1
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -29 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  43 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 43
MCLin4(1,3,1):  43 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -29 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -29
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -48 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  24 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 24
MCLin9(3,1,1):  24 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -48 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -48
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 10
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 4
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -43 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  29 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 29
MCLin4(1,3,2):  29 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -43 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -43
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -43 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  29 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 29
MCLin9(3,1,2):  29 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -43 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -43
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 1
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 1
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -37 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  32 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 32
MCLin4(2,3,0):  32 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -37 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -37
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -34 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  35 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 35
MCLin9(3,2,0):  35 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -34 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -34
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 4
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 8
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -37 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  32 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 32
MCLin4(2,3,1):  32 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -37 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -37
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -44 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  25 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 25
MCLin9(3,2,1):  25 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -44 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -44
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 10
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 3
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -39 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  30 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 30
MCLin4(2,3,2):  30 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -39 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -39
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -43 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  26 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 26
MCLin9(3,2,2):  26 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -43 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -43
makespan(0):  C - x(0,0) >= 1
makespan(1):  C - x(1,1) >= 1
makespan(2):  C - x(2,0) >= 1
makespan(3):  C - x(3,0) >= 1
cut_clique_0(0123,1):  0.25000 x(2,1) + 0.12500 x(3,1) >= 1
cut_clique_0(023,1):  0.25000 x(2,1) + 0.12500 x(3,1) >= 1
cut_clique_0(123,1):  0.25000 x(2,1) + 0.12500 x(3,1) >= 1
cut_clique_0(0123,2):  0.00880 x(0,2) + 0.02174 x(1,2) + 0.01630 x(2,2) + 0.05435 x(3,2) >= 1
cut_clique_0(023,2):  0.01195 x(0,2) + 0.03906 x(2,2) + 0.06250 x(3,2) >= 1
cut_clique_0(123,2):  0.05738 x(2,2) + 0.04918 x(3,2) >= 1
cut_clique_0(013,2):  0.01254 x(0,2) + 0.03226 x(1,2) + 0.08065 x(3,2) >= 1
cut_clique_0:  0.01559 x(0,1) + 0.01336 x(2,1) + 0.00891 x(3,1) + 0.01651 x(1,2) + 0.01238 x(2,2) + 0.04127 x(3,2) >= 1
Bounds
 0 <= C <= 48
 11 <= x(0,0) <= 47
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -32 <= v(0,1,0) <= 43
 -48 <= v(1,0,0) <= 27
 -480 <= u(0,1,0)
 -480 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -37 <= v(0,2,0) <= 35
 -37 <= v(2,0,0) <= 35
 -480 <= u(0,2,0)
 -480 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -37 <= v(0,3,0) <= 32
 -34 <= v(3,0,0) <= 35
 -480 <= u(0,3,0)
 -480 <= u(3,0,0)
 0 <= x(1,0) <= 39
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -48 <= v(1,2,0) <= 27
 -32 <= v(2,1,0) <= 43
 -480 <= u(1,2,0)
 -480 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -48 <= v(1,3,0) <= 24
 -29 <= v(3,1,0) <= 43
 -480 <= u(1,3,0)
 -480 <= u(3,1,0)
 11 <= x(2,0) <= 47
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -37 <= v(2,3,0) <= 32
 -34 <= v(3,2,0) <= 35
 -480 <= u(2,3,0)
 -480 <= u(3,2,0)
 14 <= x(3,0) <= 47
 6 <= x(0,1) <= 42
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -42 <= v(0,1,1) <= 33
 -39 <= v(1,0,1) <= 36
 -480 <= u(0,1,1)
 -480 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -38 <= v(0,2,1) <= 34
 -47 <= v(2,0,1) <= 25
 -480 <= u(0,2,1)
 -480 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -31 <= v(0,3,1) <= 38
 -47 <= v(3,0,1) <= 22
 -480 <= u(0,3,1)
 -480 <= u(3,0,1)
 8 <= x(1,1) <= 47
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -36 <= v(1,2,1) <= 39
 -48 <= v(2,1,1) <= 27
 -480 <= u(1,2,1)
 -480 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -29 <= v(1,3,1) <= 43
 -48 <= v(3,1,1) <= 24
 -480 <= u(1,3,1)
 -480 <= u(3,1,1)
 0 <= x(2,1) <= 36
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -37 <= v(2,3,1) <= 32
 -44 <= v(3,2,1) <= 25
 -480 <= u(2,3,1)
 -480 <= u(3,2,1)
 0 <= x(3,1) <= 33
 0 <= x(0,2) <= 36
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -47 <= v(0,1,2) <= 28
 -38 <= v(1,0,2) <= 37
 -480 <= u(0,1,2)
 -480 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -47 <= v(0,2,2) <= 25
 -34 <= v(2,0,2) <= 38
 -480 <= u(0,2,2)
 -480 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -47 <= v(0,3,2) <= 22
 -38 <= v(3,0,2) <= 31
 -480 <= u(0,3,2)
 -480 <= u(3,0,2)
 4 <= x(1,2) <= 43
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -43 <= v(1,2,2) <= 32
 -39 <= v(2,1,2) <= 36
 -480 <= u(1,2,2)
 -480 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -43 <= v(1,3,2) <= 29
 -43 <= v(3,1,2) <= 29
 -480 <= u(1,3,2)
 -480 <= u(3,1,2)
 8 <= x(2,2) <= 44
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -39 <= v(2,3,2) <= 30
 -43 <= v(3,2,2) <= 26
 -480 <= u(2,3,2)
 -480 <= u(3,2,2)
 4 <= x(3,2) <= 37
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
