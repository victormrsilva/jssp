\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  x(0,0) - x(0,1) >= 4
ord(0,2):  - x(0,0) + x(0,2) >= 3
ord(1,1):  x(1,1) - x(1,2) >= 5
ord(1,2):  x(1,0) - x(1,1) >= 4
ord(2,1):  x(2,0) - x(2,1) >= 10
ord(2,2):  - x(2,0) + x(2,2) >= 9
ord(3,1):  x(3,0) - x(3,1) >= 2
ord(3,2):  - x(3,0) + x(3,2) >= 2
ord(4,1):  - x(4,0) + x(4,1) >= 7
ord(4,2):  - x(4,1) + x(4,2) >= 6
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
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 5
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 3
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -84 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  63 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 63
MCLin4(0,1,0):  63 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -84 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -84
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -71 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  76 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 76
MCLin9(1,0,0):  76 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -71 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -71
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 4
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 4
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -83 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  64 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 64
MCLin4(0,1,1):  64 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -83 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -83
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -72 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  75 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 75
MCLin9(1,0,1):  75 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -72 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -72
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 5
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 8
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -72 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  75 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 75
MCLin4(0,1,2):  75 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -72 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -72
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -88 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  59 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 59
MCLin9(1,0,2):  59 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -88 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -88
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 9
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 3
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -75 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  58 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 58
MCLin4(0,2,0):  58 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -75 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -75
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -70 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  63 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 63
MCLin9(2,0,0):  63 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -70 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -70
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 10
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 4
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -70 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  63 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 63
MCLin4(0,2,1):  63 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -70 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -70
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -77 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  56 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 56
MCLin9(2,0,1):  56 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -77 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -77
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 9
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 8
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -81 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  52 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 52
MCLin4(0,2,2):  52 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -81 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -81
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -69 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  64 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 64
MCLin9(2,0,2):  64 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -69 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -69
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 2
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 3
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -74 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  73 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 73
MCLin4(0,3,0):  73 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -74 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -74
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -78 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  69 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 69
MCLin9(3,0,0):  69 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -78 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -78
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 2
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 4
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -76 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  71 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 71
MCLin4(0,3,1):  71 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -76 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -76
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -77 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  70 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 70
MCLin9(3,0,1):  70 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -77 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -77
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 10
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 8
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -81 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  66 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 66
MCLin4(0,3,2):  66 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -81 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -81
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -84 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  63 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 63
MCLin9(3,0,2):  63 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -84 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -84
MCLinV(0,4,0):  x(0,0) - v(0,4,0) - x(4,0) = 7
MCLinV(4,0,0):  - x(0,0) - v(4,0,0) + x(4,0) = 3
MCLin1(0,4,0):  u(0,4,0) >= -0
MCLin2(0,4,0):  -74 y(4,0,0) - u(0,4,0) <= -0
MCLin3(0,4,0):  70 y(4,0,0) + v(0,4,0) - u(0,4,0) <= 70
MCLin4(0,4,0):  70 y(4,0,0) - u(0,4,0) >= -0
MCLin5(0,4,0):  -74 y(4,0,0) + v(0,4,0) - u(0,4,0) >= -74
MCLin6(4,0,0):  u(4,0,0) >= -0
MCLin7(4,0,0):  -80 y(0,4,0) - u(4,0,0) <= -0
MCLin8(4,0,0):  64 y(0,4,0) + v(4,0,0) - u(4,0,0) <= 64
MCLin9(4,0,0):  64 y(0,4,0) - u(4,0,0) >= -0
MCLin10(4,0,0):  -80 y(0,4,0) + v(4,0,0) - u(4,0,0) >= -80
MCLinV(0,4,1):  x(0,1) - v(0,4,1) - x(4,1) = 6
MCLinV(4,0,1):  - x(0,1) - v(4,0,1) + x(4,1) = 4
MCLin1(0,4,1):  u(0,4,1) >= -0
MCLin2(0,4,1):  -84 y(4,0,1) - u(0,4,1) <= -0
MCLin3(0,4,1):  60 y(4,0,1) + v(0,4,1) - u(0,4,1) <= 60
MCLin4(0,4,1):  60 y(4,0,1) - u(0,4,1) >= -0
MCLin5(0,4,1):  -84 y(4,0,1) + v(0,4,1) - u(0,4,1) >= -84
MCLin6(4,0,1):  u(4,0,1) >= -0
MCLin7(4,0,1):  -70 y(0,4,1) - u(4,0,1) <= -0
MCLin8(4,0,1):  74 y(0,4,1) + v(4,0,1) - u(4,0,1) <= 74
MCLin9(4,0,1):  74 y(0,4,1) - u(4,0,1) >= -0
MCLin10(4,0,1):  -70 y(0,4,1) + v(4,0,1) - u(4,0,1) >= -70
MCLinV(0,4,2):  x(0,2) - v(0,4,2) - x(4,2) = 4
MCLinV(4,0,2):  - x(0,2) - v(4,0,2) + x(4,2) = 8
MCLin1(0,4,2):  u(0,4,2) >= -0
MCLin2(0,4,2):  -81 y(4,0,2) - u(0,4,2) <= -0
MCLin3(0,4,2):  63 y(4,0,2) + v(0,4,2) - u(0,4,2) <= 63
MCLin4(0,4,2):  63 y(4,0,2) - u(0,4,2) >= -0
MCLin5(0,4,2):  -81 y(4,0,2) + v(0,4,2) - u(0,4,2) >= -81
MCLin6(4,0,2):  u(4,0,2) >= -0
MCLin7(4,0,2):  -75 y(0,4,2) - u(4,0,2) <= -0
MCLin8(4,0,2):  69 y(0,4,2) + v(4,0,2) - u(4,0,2) <= 69
MCLin9(4,0,2):  69 y(0,4,2) - u(4,0,2) >= -0
MCLin10(4,0,2):  -75 y(0,4,2) + v(4,0,2) - u(4,0,2) >= -75
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 9
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 5
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -70 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  64 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 64
MCLin4(1,2,0):  64 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -70 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -70
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -78 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  56 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 56
MCLin9(2,1,0):  56 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -78 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -78
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 10
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 4
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -65 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  69 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 69
MCLin4(1,2,1):  69 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -65 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -65
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -83 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  51 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 51
MCLin9(2,1,1):  51 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -83 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -83
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 9
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 5
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -88 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  46 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 46
MCLin4(1,2,2):  46 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -88 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -88
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -60 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  74 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 74
MCLin9(2,1,2):  74 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -60 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -60
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 2
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 5
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -69 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  79 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 79
MCLin4(1,3,0):  79 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -69 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -69
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -86 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  62 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 62
MCLin9(3,1,0):  62 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -86 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -86
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 2
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 4
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -71 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  77 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 77
MCLin4(1,3,1):  77 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -71 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -71
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -83 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  65 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 65
MCLin9(3,1,1):  65 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -83 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -83
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 10
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 5
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -88 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  60 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 60
MCLin4(1,3,2):  60 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -88 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -88
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -75 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  73 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 73
MCLin9(3,1,2):  73 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -75 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -75
MCLinV(1,4,0):  x(1,0) - v(1,4,0) - x(4,0) = 7
MCLinV(4,1,0):  - x(1,0) - v(4,1,0) + x(4,0) = 5
MCLin1(1,4,0):  u(1,4,0) >= -0
MCLin2(1,4,0):  -69 y(4,1,0) - u(1,4,0) <= -0
MCLin3(1,4,0):  76 y(4,1,0) + v(1,4,0) - u(1,4,0) <= 76
MCLin4(1,4,0):  76 y(4,1,0) - u(1,4,0) >= -0
MCLin5(1,4,0):  -69 y(4,1,0) + v(1,4,0) - u(1,4,0) >= -69
MCLin6(4,1,0):  u(4,1,0) >= -0
MCLin7(4,1,0):  -88 y(1,4,0) - u(4,1,0) <= -0
MCLin8(4,1,0):  57 y(1,4,0) + v(4,1,0) - u(4,1,0) <= 57
MCLin9(4,1,0):  57 y(1,4,0) - u(4,1,0) >= -0
MCLin10(4,1,0):  -88 y(1,4,0) + v(4,1,0) - u(4,1,0) >= -88
MCLinV(1,4,1):  x(1,1) - v(1,4,1) - x(4,1) = 6
MCLinV(4,1,1):  - x(1,1) - v(4,1,1) + x(4,1) = 4
MCLin1(1,4,1):  u(1,4,1) >= -0
MCLin2(1,4,1):  -79 y(4,1,1) - u(1,4,1) <= -0
MCLin3(1,4,1):  66 y(4,1,1) + v(1,4,1) - u(1,4,1) <= 66
MCLin4(1,4,1):  66 y(4,1,1) - u(1,4,1) >= -0
MCLin5(1,4,1):  -79 y(4,1,1) + v(1,4,1) - u(1,4,1) >= -79
MCLin6(4,1,1):  u(4,1,1) >= -0
MCLin7(4,1,1):  -76 y(1,4,1) - u(4,1,1) <= -0
MCLin8(4,1,1):  69 y(1,4,1) + v(4,1,1) - u(4,1,1) <= 69
MCLin9(4,1,1):  69 y(1,4,1) - u(4,1,1) >= -0
MCLin10(4,1,1):  -76 y(1,4,1) + v(4,1,1) - u(4,1,1) >= -76
MCLinV(1,4,2):  x(1,2) - v(1,4,2) - x(4,2) = 4
MCLinV(4,1,2):  - x(1,2) - v(4,1,2) + x(4,2) = 5
MCLin1(1,4,2):  u(1,4,2) >= -0
MCLin2(1,4,2):  -88 y(4,1,2) - u(1,4,2) <= -0
MCLin3(1,4,2):  57 y(4,1,2) + v(1,4,2) - u(1,4,2) <= 57
MCLin4(1,4,2):  57 y(4,1,2) - u(1,4,2) >= -0
MCLin5(1,4,2):  -88 y(4,1,2) + v(1,4,2) - u(1,4,2) >= -88
MCLin6(4,1,2):  u(4,1,2) >= -0
MCLin7(4,1,2):  -66 y(1,4,2) - u(4,1,2) <= -0
MCLin8(4,1,2):  79 y(1,4,2) + v(4,1,2) - u(4,1,2) <= 79
MCLin9(4,1,2):  79 y(1,4,2) - u(4,1,2) >= -0
MCLin10(4,1,2):  -66 y(1,4,2) + v(4,1,2) - u(4,1,2) >= -66
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 2
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 9
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -68 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  66 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 66
MCLin4(2,3,0):  66 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -68 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -68
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -77 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  57 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 57
MCLin9(3,2,0):  57 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -77 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -77
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 2
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 10
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -76 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  58 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 58
MCLin4(2,3,1):  58 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -76 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -76
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -70 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  64 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 64
MCLin9(3,2,1):  64 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -70 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -70
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 10
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 9
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -69 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  65 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 65
MCLin4(2,3,2):  65 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -69 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -69
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -84 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  50 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 50
MCLin9(3,2,2):  50 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -84 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -84
MCLinV(2,4,0):  x(2,0) - v(2,4,0) - x(4,0) = 7
MCLinV(4,2,0):  - x(2,0) - v(4,2,0) + x(4,0) = 9
MCLin1(2,4,0):  u(2,4,0) >= -0
MCLin2(2,4,0):  -68 y(4,2,0) - u(2,4,0) <= -0
MCLin3(2,4,0):  63 y(4,2,0) + v(2,4,0) - u(2,4,0) <= 63
MCLin4(2,4,0):  63 y(4,2,0) - u(2,4,0) >= -0
MCLin5(2,4,0):  -68 y(4,2,0) + v(2,4,0) - u(2,4,0) >= -68
MCLin6(4,2,0):  u(4,2,0) >= -0
MCLin7(4,2,0):  -79 y(2,4,0) - u(4,2,0) <= -0
MCLin8(4,2,0):  52 y(2,4,0) + v(4,2,0) - u(4,2,0) <= 52
MCLin9(4,2,0):  52 y(2,4,0) - u(4,2,0) >= -0
MCLin10(4,2,0):  -79 y(2,4,0) + v(4,2,0) - u(4,2,0) >= -79
MCLinV(2,4,1):  x(2,1) - v(2,4,1) - x(4,1) = 6
MCLinV(4,2,1):  - x(2,1) - v(4,2,1) + x(4,1) = 10
MCLin1(2,4,1):  u(2,4,1) >= -0
MCLin2(2,4,1):  -84 y(4,2,1) - u(2,4,1) <= -0
MCLin3(2,4,1):  47 y(4,2,1) + v(2,4,1) - u(2,4,1) <= 47
MCLin4(2,4,1):  47 y(4,2,1) - u(2,4,1) >= -0
MCLin5(2,4,1):  -84 y(4,2,1) + v(2,4,1) - u(2,4,1) >= -84
MCLin6(4,2,1):  u(4,2,1) >= -0
MCLin7(4,2,1):  -63 y(2,4,1) - u(4,2,1) <= -0
MCLin8(4,2,1):  68 y(2,4,1) + v(4,2,1) - u(4,2,1) <= 68
MCLin9(4,2,1):  68 y(2,4,1) - u(4,2,1) >= -0
MCLin10(4,2,1):  -63 y(2,4,1) + v(4,2,1) - u(4,2,1) >= -63
MCLinV(2,4,2):  x(2,2) - v(2,4,2) - x(4,2) = 4
MCLinV(4,2,2):  - x(2,2) - v(4,2,2) + x(4,2) = 9
MCLin1(2,4,2):  u(2,4,2) >= -0
MCLin2(2,4,2):  -69 y(4,2,2) - u(2,4,2) <= -0
MCLin3(2,4,2):  62 y(4,2,2) + v(2,4,2) - u(2,4,2) <= 62
MCLin4(2,4,2):  62 y(4,2,2) - u(2,4,2) >= -0
MCLin5(2,4,2):  -69 y(4,2,2) + v(2,4,2) - u(2,4,2) >= -69
MCLin6(4,2,2):  u(4,2,2) >= -0
MCLin7(4,2,2):  -75 y(2,4,2) - u(4,2,2) <= -0
MCLin8(4,2,2):  56 y(2,4,2) + v(4,2,2) - u(4,2,2) <= 56
MCLin9(4,2,2):  56 y(2,4,2) - u(4,2,2) >= -0
MCLin10(4,2,2):  -75 y(2,4,2) + v(4,2,2) - u(4,2,2) >= -75
MCLinV(3,4,0):  x(3,0) - v(3,4,0) - x(4,0) = 7
MCLinV(4,3,0):  - x(3,0) - v(4,3,0) + x(4,0) = 2
MCLin1(3,4,0):  u(3,4,0) >= -0
MCLin2(3,4,0):  -76 y(4,3,0) - u(3,4,0) <= -0
MCLin3(3,4,0):  69 y(4,3,0) + v(3,4,0) - u(3,4,0) <= 69
MCLin4(3,4,0):  69 y(4,3,0) - u(3,4,0) >= -0
MCLin5(3,4,0):  -76 y(4,3,0) + v(3,4,0) - u(3,4,0) >= -76
MCLin6(4,3,0):  u(4,3,0) >= -0
MCLin7(4,3,0):  -78 y(3,4,0) - u(4,3,0) <= -0
MCLin8(4,3,0):  67 y(3,4,0) + v(4,3,0) - u(4,3,0) <= 67
MCLin9(4,3,0):  67 y(3,4,0) - u(4,3,0) >= -0
MCLin10(4,3,0):  -78 y(3,4,0) + v(4,3,0) - u(4,3,0) >= -78
MCLinV(3,4,1):  x(3,1) - v(3,4,1) - x(4,1) = 6
MCLinV(4,3,1):  - x(3,1) - v(4,3,1) + x(4,1) = 2
MCLin1(3,4,1):  u(3,4,1) >= -0
MCLin2(3,4,1):  -84 y(4,3,1) - u(3,4,1) <= -0
MCLin3(3,4,1):  61 y(4,3,1) + v(3,4,1) - u(3,4,1) <= 61
MCLin4(3,4,1):  61 y(4,3,1) - u(3,4,1) >= -0
MCLin5(3,4,1):  -84 y(4,3,1) + v(3,4,1) - u(3,4,1) >= -84
MCLin6(4,3,1):  u(4,3,1) >= -0
MCLin7(4,3,1):  -69 y(3,4,1) - u(4,3,1) <= -0
MCLin8(4,3,1):  76 y(3,4,1) + v(4,3,1) - u(4,3,1) <= 76
MCLin9(4,3,1):  76 y(3,4,1) - u(4,3,1) >= -0
MCLin10(4,3,1):  -69 y(3,4,1) + v(4,3,1) - u(4,3,1) >= -69
MCLinV(3,4,2):  x(3,2) - v(3,4,2) - x(4,2) = 4
MCLinV(4,3,2):  - x(3,2) - v(4,3,2) + x(4,2) = 10
MCLin1(3,4,2):  u(3,4,2) >= -0
MCLin2(3,4,2):  -84 y(4,3,2) - u(3,4,2) <= -0
MCLin3(3,4,2):  61 y(4,3,2) + v(3,4,2) - u(3,4,2) <= 61
MCLin4(3,4,2):  61 y(4,3,2) - u(3,4,2) >= -0
MCLin5(3,4,2):  -84 y(4,3,2) + v(3,4,2) - u(3,4,2) >= -84
MCLin6(4,3,2):  u(4,3,2) >= -0
MCLin7(4,3,2):  -75 y(3,4,2) - u(4,3,2) <= -0
MCLin8(4,3,2):  70 y(3,4,2) + v(4,3,2) - u(4,3,2) <= 70
MCLin9(4,3,2):  70 y(3,4,2) - u(4,3,2) >= -0
MCLin10(4,3,2):  -75 y(3,4,2) + v(4,3,2) - u(4,3,2) >= -75
makespan(0):  C - x(0,2) >= 8
makespan(1):  C - x(1,0) >= 5
makespan(2):  C - x(2,2) >= 9
makespan(3):  C - x(3,2) >= 10
makespan(4):  C - x(4,2) >= 4
cut_clique_0(01234,0):  0.01181 x(0,0) + 0.01969 x(1,0) + 0.03543 x(2,0) + 0.00787 x(3,0) + 0.01147 x(4,0) >= 1
cut_clique_0(0123,0):  0.02273 x(1,0) + 0.05682 x(2,0) >= 1
cut_clique_0(0234,0):  0.01316 x(0,0) + 0.06579 x(2,0) + 0.01316 x(3,0) + 0.01385 x(4,0) >= 1
cut_clique_0(0124,0):  0.00458 x(0,0) + 0.02616 x(1,0) + 0.04709 x(2,0) + 0.00378 x(4,0) >= 1
cut_clique_0(1234,0):  0.02273 x(1,0) + 0.05682 x(2,0) >= 1
cut_clique_0(012,0):  0.02273 x(1,0) + 0.05682 x(2,0) >= 1
cut_clique_0(123,0):  0.02273 x(1,0) + 0.05682 x(2,0) >= 1
cut_clique_0(124,0):  0.02273 x(1,0) + 0.05682 x(2,0) >= 1
cut_clique_0(01234,1):  0.02439 x(0,1) + 0.06098 x(2,1) + 0.01220 x(3,1) + 0.03659 x(4,1) >= 1
cut_clique_0(0134,1):  0.02439 x(1,1) + 0.09756 x(4,1) >= 1
cut_clique_0(0234,1):  0.02439 x(0,1) + 0.06098 x(2,1) + 0.01220 x(3,1) + 0.03659 x(4,1) >= 1
cut_clique_0(0123,1):  0.05882 x(0,1) + 0.14706 x(2,1) + 0.02941 x(3,1) >= 1
cut_clique_0(0124,1):  0.10000 x(0,1) + 0.25000 x(2,1) >= 1
cut_clique_0(1234,1):  0.50000 x(2,1) + 0.10000 x(3,1) >= 1
cut_clique_0(134,1):  0.02439 x(1,1) + 0.09756 x(4,1) >= 1
cut_clique_0(014,1):  0.02439 x(1,1) + 0.09756 x(4,1) >= 1
cut_clique_0(023,1):  0.05882 x(0,1) + 0.14706 x(2,1) + 0.02941 x(3,1) >= 1
cut_clique_0(024,1):  0.10000 x(0,1) + 0.25000 x(2,1) >= 1
cut_clique_0(234,1):  0.50000 x(2,1) + 0.10000 x(3,1) >= 1
cut_clique_0(012,1):  0.10000 x(0,1) + 0.25000 x(2,1) >= 1
cut_clique_0(123,1):  0.50000 x(2,1) + 0.10000 x(3,1) >= 1
cut_clique_0(124,1):  0.03226 x(1,1) + 0.02688 x(2,1) + 0.04839 x(4,1) >= 1
cut_clique_0(01234,2):  0.00453 x(0,2) + 0.00183 x(1,2) + 0.00510 x(2,2) + 0.00414 x(3,2) + 0.05037 x(4,2) >= 1
cut_clique_0(0234,2):  0.00488 x(0,2) + 0.00697 x(2,2) + 0.00212 x(3,2) + 0.05162 x(4,2) >= 1
cut_clique_0(0134,2):  0.00506 x(0,2) + 0.00267 x(1,2) + 0.00632 x(3,2) + 0.05815 x(4,2) >= 1
cut_clique_0(1234,2):  0.00219 x(1,2) + 0.00675 x(2,2) + 0.00491 x(3,2) + 0.05648 x(4,2) >= 1
cut_clique_0(0124,2):  0.00502 x(0,2) + 0.00920 x(2,2) + 0.05268 x(4,2) >= 1
cut_clique_0(0123,2):  0.01060 x(0,2) + 0.00177 x(1,2) + 0.03446 x(2,2) + 0.00968 x(3,2) >= 1
cut_clique_0(024,2):  0.00502 x(0,2) + 0.00920 x(2,2) + 0.05268 x(4,2) >= 1
cut_clique_0(134,2):  0.00409 x(1,2) + 0.01053 x(3,2) + 0.06316 x(4,2) >= 1
cut_clique_0(034,2):  0.00266 x(0,2) + 0.00332 x(3,2) + 0.06707 x(4,2) >= 1
cut_clique_0(023,2):  0.00832 x(0,2) + 0.03883 x(2,2) + 0.00728 x(3,2) >= 1
cut_clique_0(234,2):  0.00376 x(2,2) + 0.00273 x(3,2) + 0.06555 x(4,2) >= 1
cut_clique_1(01234,1):  0.01587 x(0,1) + 0.01587 x(1,1) + 0.03968 x(2,1) + 0.00794 x(3,1) + 0.02381 x(4,1) >= 1
cut_clique_0:  0.03243 x(2,1) + 0.04865 x(2,2) + 0.00270 x(3,2) >= 1
cut_clique_1:  0.00203 x(4,0) + 0.02838 x(2,1) + 0.04865 x(2,2) + 0.00270 x(3,2) >= 1
cut_clique_2:  0.01571 x(0,0) + 0.04712 x(2,0) + 0.01047 x(3,0) + 0.00238 x(4,0) + 0.07853 x(2,1) >= 1
Bounds
 0 <= C <= 88
 4 <= x(0,0) <= 77
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -84 <= v(0,1,0) <= 63
 -71 <= v(1,0,0) <= 76
 -880 <= u(0,1,0)
 -880 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -75 <= v(0,2,0) <= 58
 -70 <= v(2,0,0) <= 63
 -880 <= u(0,2,0)
 -880 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -74 <= v(0,3,0) <= 73
 -78 <= v(3,0,0) <= 69
 -880 <= u(0,3,0)
 -880 <= u(3,0,0)
 0 <= y(0,4,0) <= 1
 0 <= y(4,0,0) <= 1
 -74 <= v(0,4,0) <= 70
 -80 <= v(4,0,0) <= 64
 -880 <= u(0,4,0)
 -880 <= u(4,0,0)
 9 <= x(1,0) <= 83
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -70 <= v(1,2,0) <= 64
 -78 <= v(2,1,0) <= 56
 -880 <= u(1,2,0)
 -880 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -69 <= v(1,3,0) <= 79
 -86 <= v(3,1,0) <= 62
 -880 <= u(1,3,0)
 -880 <= u(3,1,0)
 0 <= y(1,4,0) <= 1
 0 <= y(4,1,0) <= 1
 -69 <= v(1,4,0) <= 76
 -88 <= v(4,1,0) <= 57
 -880 <= u(1,4,0)
 -880 <= u(4,1,0)
 10 <= x(2,0) <= 70
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -68 <= v(2,3,0) <= 66
 -77 <= v(3,2,0) <= 57
 -880 <= u(2,3,0)
 -880 <= u(3,2,0)
 0 <= y(2,4,0) <= 1
 0 <= y(4,2,0) <= 1
 -68 <= v(2,4,0) <= 63
 -79 <= v(4,2,0) <= 52
 -880 <= u(2,4,0)
 -880 <= u(4,2,0)
 2 <= x(3,0) <= 76
 0 <= y(3,4,0) <= 1
 0 <= y(4,3,0) <= 1
 -76 <= v(3,4,0) <= 69
 -78 <= v(4,3,0) <= 67
 -880 <= u(3,4,0)
 -880 <= u(4,3,0)
 0 <= x(4,0) <= 71
 0 <= x(0,1) <= 73
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -83 <= v(0,1,1) <= 64
 -72 <= v(1,0,1) <= 75
 -880 <= u(0,1,1)
 -880 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -70 <= v(0,2,1) <= 63
 -77 <= v(2,0,1) <= 56
 -880 <= u(0,2,1)
 -880 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -76 <= v(0,3,1) <= 71
 -77 <= v(3,0,1) <= 70
 -880 <= u(0,3,1)
 -880 <= u(3,0,1)
 0 <= y(0,4,1) <= 1
 0 <= y(4,0,1) <= 1
 -84 <= v(0,4,1) <= 60
 -70 <= v(4,0,1) <= 74
 -880 <= u(0,4,1)
 -880 <= u(4,0,1)
 5 <= x(1,1) <= 79
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -65 <= v(1,2,1) <= 69
 -83 <= v(2,1,1) <= 51
 -880 <= u(1,2,1)
 -880 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -71 <= v(1,3,1) <= 77
 -83 <= v(3,1,1) <= 65
 -880 <= u(1,3,1)
 -880 <= u(3,1,1)
 0 <= y(1,4,1) <= 1
 0 <= y(4,1,1) <= 1
 -79 <= v(1,4,1) <= 66
 -76 <= v(4,1,1) <= 69
 -880 <= u(1,4,1)
 -880 <= u(4,1,1)
 0 <= x(2,1) <= 60
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -76 <= v(2,3,1) <= 58
 -70 <= v(3,2,1) <= 64
 -880 <= u(2,3,1)
 -880 <= u(3,2,1)
 0 <= y(2,4,1) <= 1
 0 <= y(4,2,1) <= 1
 -84 <= v(2,4,1) <= 47
 -63 <= v(4,2,1) <= 68
 -880 <= u(2,4,1)
 -880 <= u(4,2,1)
 0 <= x(3,1) <= 74
 0 <= y(3,4,1) <= 1
 0 <= y(4,3,1) <= 1
 -84 <= v(3,4,1) <= 61
 -69 <= v(4,3,1) <= 76
 -880 <= u(3,4,1)
 -880 <= u(4,3,1)
 7 <= x(4,1) <= 78
 7 <= x(0,2) <= 80
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -72 <= v(0,1,2) <= 75
 -88 <= v(1,0,2) <= 59
 -880 <= u(0,1,2)
 -880 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -81 <= v(0,2,2) <= 52
 -69 <= v(2,0,2) <= 64
 -880 <= u(0,2,2)
 -880 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -81 <= v(0,3,2) <= 66
 -84 <= v(3,0,2) <= 63
 -880 <= u(0,3,2)
 -880 <= u(3,0,2)
 0 <= y(0,4,2) <= 1
 0 <= y(4,0,2) <= 1
 -81 <= v(0,4,2) <= 63
 -75 <= v(4,0,2) <= 69
 -880 <= u(0,4,2)
 -880 <= u(4,0,2)
 0 <= x(1,2) <= 74
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -88 <= v(1,2,2) <= 46
 -60 <= v(2,1,2) <= 74
 -880 <= u(1,2,2)
 -880 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -88 <= v(1,3,2) <= 60
 -75 <= v(3,1,2) <= 73
 -880 <= u(1,3,2)
 -880 <= u(3,1,2)
 0 <= y(1,4,2) <= 1
 0 <= y(4,1,2) <= 1
 -88 <= v(1,4,2) <= 57
 -66 <= v(4,1,2) <= 79
 -880 <= u(1,4,2)
 -880 <= u(4,1,2)
 19 <= x(2,2) <= 79
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -69 <= v(2,3,2) <= 65
 -84 <= v(3,2,2) <= 50
 -880 <= u(2,3,2)
 -880 <= u(3,2,2)
 0 <= y(2,4,2) <= 1
 0 <= y(4,2,2) <= 1
 -69 <= v(2,4,2) <= 62
 -75 <= v(4,2,2) <= 56
 -880 <= u(2,4,2)
 -880 <= u(4,2,2)
 4 <= x(3,2) <= 78
 0 <= y(3,4,2) <= 1
 0 <= y(4,3,2) <= 1
 -84 <= v(3,4,2) <= 61
 -75 <= v(4,3,2) <= 70
 -880 <= u(3,4,2)
 -880 <= u(4,3,2)
 13 <= x(4,2) <= 84
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
