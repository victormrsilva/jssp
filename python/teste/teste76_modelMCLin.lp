\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,2) >= 5
ord(0,2):  x(0,0) - x(0,2) >= 5
ord(0,3):  - x(0,0) + x(0,3) >= 8
ord(1,1):  - x(1,2) + x(1,3) >= 2
ord(1,2):  x(1,0) - x(1,3) >= 2
ord(1,3):  - x(1,0) + x(1,1) >= 4
ord(2,1):  - x(2,1) + x(2,2) >= 8
ord(2,2):  x(2,0) - x(2,2) >= 10
ord(2,3):  - x(2,0) + x(2,3) >= 3
ord(3,1):  x(3,1) - x(3,2) >= 9
ord(3,2):  - x(3,1) + x(3,3) >= 6
ord(3,3):  x(3,0) - x(3,3) >= 7
ord(4,1):  x(4,0) - x(4,2) >= 2
ord(4,2):  - x(4,0) + x(4,3) >= 1
ord(4,3):  x(4,1) - x(4,3) >= 4
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,1,3):  y(0,1,3) + y(1,0,3) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(0,2,3):  y(0,2,3) + y(2,0,3) = 1
triangle2(0,3,0):  y(0,3,0) + y(3,0,0) = 1
triangle2(0,3,1):  y(0,3,1) + y(3,0,1) = 1
triangle2(0,3,2):  y(0,3,2) + y(3,0,2) = 1
triangle2(0,3,3):  y(0,3,3) + y(3,0,3) = 1
triangle2(0,4,0):  y(0,4,0) + y(4,0,0) = 1
triangle2(0,4,1):  y(0,4,1) + y(4,0,1) = 1
triangle2(0,4,2):  y(0,4,2) + y(4,0,2) = 1
triangle2(0,4,3):  y(0,4,3) + y(4,0,3) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
triangle2(1,2,3):  y(1,2,3) + y(2,1,3) = 1
triangle2(1,3,0):  y(1,3,0) + y(3,1,0) = 1
triangle2(1,3,1):  y(1,3,1) + y(3,1,1) = 1
triangle2(1,3,2):  y(1,3,2) + y(3,1,2) = 1
triangle2(1,3,3):  y(1,3,3) + y(3,1,3) = 1
triangle2(1,4,0):  y(1,4,0) + y(4,1,0) = 1
triangle2(1,4,1):  y(1,4,1) + y(4,1,1) = 1
triangle2(1,4,2):  y(1,4,2) + y(4,1,2) = 1
triangle2(1,4,3):  y(1,4,3) + y(4,1,3) = 1
triangle2(2,3,0):  y(2,3,0) + y(3,2,0) = 1
triangle2(2,3,1):  y(2,3,1) + y(3,2,1) = 1
triangle2(2,3,2):  y(2,3,2) + y(3,2,2) = 1
triangle2(2,3,3):  y(2,3,3) + y(3,2,3) = 1
triangle2(2,4,0):  y(2,4,0) + y(4,2,0) = 1
triangle2(2,4,1):  y(2,4,1) + y(4,2,1) = 1
triangle2(2,4,2):  y(2,4,2) + y(4,2,2) = 1
triangle2(2,4,3):  y(2,4,3) + y(4,2,3) = 1
triangle2(3,4,0):  y(3,4,0) + y(4,3,0) = 1
triangle2(3,4,1):  y(3,4,1) + y(4,3,1) = 1
triangle2(3,4,2):  y(3,4,2) + y(4,3,2) = 1
triangle2(3,4,3):  y(3,4,3) + y(4,3,3) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 4
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 8
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -83 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  72 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 72
MCLin4(0,1,0):  72 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -83 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -83
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -84 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  71 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 71
MCLin9(1,0,0):  71 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -84 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -84
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 2
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 5
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -95 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  60 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 60
MCLin4(0,1,1):  60 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -95 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -95
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -67 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  88 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 88
MCLin9(1,0,1):  88 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -67 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -67
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 2
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 5
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -82 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  73 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 73
MCLin4(0,1,2):  73 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -82 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -82
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -80 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  75 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 75
MCLin9(1,0,2):  75 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -80 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -80
MCLinV(0,1,3):  x(0,3) - v(0,1,3) - x(1,3) = 2
MCLinV(1,0,3):  - x(0,3) - v(1,0,3) + x(1,3) = 7
MCLin1(0,1,3):  u(0,1,3) >= -0
MCLin2(0,1,3):  -71 y(1,0,3) - u(0,1,3) <= -0
MCLin3(0,1,3):  84 y(1,0,3) + v(0,1,3) - u(0,1,3) <= 84
MCLin4(0,1,3):  84 y(1,0,3) - u(0,1,3) >= -0
MCLin5(0,1,3):  -71 y(1,0,3) + v(0,1,3) - u(0,1,3) >= -71
MCLin6(1,0,3):  u(1,0,3) >= -0
MCLin7(1,0,3):  -93 y(0,1,3) - u(1,0,3) <= -0
MCLin8(1,0,3):  62 y(0,1,3) + v(1,0,3) - u(1,0,3) <= 62
MCLin9(1,0,3):  62 y(0,1,3) - u(1,0,3) >= -0
MCLin10(1,0,3):  -93 y(0,1,3) + v(1,0,3) - u(1,0,3) >= -93
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 3
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 8
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -83 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  59 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 59
MCLin4(0,2,0):  59 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -83 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -83
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -70 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  72 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 72
MCLin9(2,0,0):  72 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -70 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -70
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 8
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 5
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -80 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  62 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 62
MCLin4(0,2,1):  62 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -80 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -80
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -75 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  67 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 67
MCLin9(2,0,1):  67 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -75 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -75
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 10
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 5
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -85 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  57 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 57
MCLin4(0,2,2):  57 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -85 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -85
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -72 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  70 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 70
MCLin9(2,0,2):  70 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -72 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -72
MCLinV(0,2,3):  x(0,3) - v(0,2,3) - x(2,3) = 2
MCLinV(2,0,3):  - x(0,3) - v(2,0,3) + x(2,3) = 7
MCLin1(0,2,3):  u(0,2,3) >= -0
MCLin2(0,2,3):  -77 y(2,0,3) - u(0,2,3) <= -0
MCLin3(0,2,3):  65 y(2,0,3) + v(0,2,3) - u(0,2,3) <= 65
MCLin4(0,2,3):  65 y(2,0,3) - u(0,2,3) >= -0
MCLin5(0,2,3):  -77 y(2,0,3) + v(0,2,3) - u(0,2,3) >= -77
MCLin6(2,0,3):  u(2,0,3) >= -0
MCLin7(2,0,3):  -74 y(0,2,3) - u(2,0,3) <= -0
MCLin8(2,0,3):  68 y(0,2,3) + v(2,0,3) - u(2,0,3) <= 68
MCLin9(2,0,3):  68 y(0,2,3) - u(2,0,3) >= -0
MCLin10(2,0,3):  -74 y(0,2,3) + v(2,0,3) - u(2,0,3) >= -74
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 5
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 8
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -85 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  53 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 53
MCLin4(0,3,0):  53 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -85 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -85
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -66 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  72 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 72
MCLin9(3,0,0):  72 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -66 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -66
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 6
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 5
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -83 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  55 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 55
MCLin4(0,3,1):  55 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -83 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -83
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -66 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  72 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 72
MCLin9(3,0,1):  72 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -66 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -66
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 9
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 5
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -72 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  66 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 66
MCLin4(0,3,2):  66 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -72 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -72
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -80 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  58 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 58
MCLin9(3,0,2):  58 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -80 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -80
MCLinV(0,3,3):  x(0,3) - v(0,3,3) - x(3,3) = 7
MCLinV(3,0,3):  - x(0,3) - v(3,0,3) + x(3,3) = 7
MCLin1(0,3,3):  u(0,3,3) >= -0
MCLin2(0,3,3):  -72 y(3,0,3) - u(0,3,3) <= -0
MCLin3(0,3,3):  66 y(3,0,3) + v(0,3,3) - u(0,3,3) <= 66
MCLin4(0,3,3):  66 y(3,0,3) - u(0,3,3) >= -0
MCLin5(0,3,3):  -72 y(3,0,3) + v(0,3,3) - u(0,3,3) >= -72
MCLin6(3,0,3):  u(3,0,3) >= -0
MCLin7(3,0,3):  -80 y(0,3,3) - u(3,0,3) <= -0
MCLin8(3,0,3):  58 y(0,3,3) + v(3,0,3) - u(3,0,3) <= 58
MCLin9(3,0,3):  58 y(0,3,3) - u(3,0,3) >= -0
MCLin10(3,0,3):  -80 y(0,3,3) + v(3,0,3) - u(3,0,3) >= -80
MCLinV(0,4,0):  x(0,0) - v(0,4,0) - x(4,0) = 1
MCLinV(4,0,0):  - x(0,0) - v(4,0,0) + x(4,0) = 8
MCLin1(0,4,0):  u(0,4,0) >= -0
MCLin2(0,4,0):  -78 y(4,0,0) - u(0,4,0) <= -0
MCLin3(0,4,0):  77 y(4,0,0) + v(0,4,0) - u(0,4,0) <= 77
MCLin4(0,4,0):  77 y(4,0,0) - u(0,4,0) >= -0
MCLin5(0,4,0):  -78 y(4,0,0) + v(0,4,0) - u(0,4,0) >= -78
MCLin6(4,0,0):  u(4,0,0) >= -0
MCLin7(4,0,0):  -86 y(0,4,0) - u(4,0,0) <= -0
MCLin8(4,0,0):  69 y(0,4,0) + v(4,0,0) - u(4,0,0) <= 69
MCLin9(4,0,0):  69 y(0,4,0) - u(4,0,0) >= -0
MCLin10(4,0,0):  -86 y(0,4,0) + v(4,0,0) - u(4,0,0) >= -86
MCLinV(0,4,1):  x(0,1) - v(0,4,1) - x(4,1) = 3
MCLinV(4,0,1):  - x(0,1) - v(4,0,1) + x(4,1) = 5
MCLin1(0,4,1):  u(0,4,1) >= -0
MCLin2(0,4,1):  -95 y(4,0,1) - u(0,4,1) <= -0
MCLin3(0,4,1):  60 y(4,0,1) + v(0,4,1) - u(0,4,1) <= 60
MCLin4(0,4,1):  60 y(4,0,1) - u(0,4,1) >= -0
MCLin5(0,4,1):  -95 y(4,0,1) + v(0,4,1) - u(0,4,1) >= -95
MCLin6(4,0,1):  u(4,0,1) >= -0
MCLin7(4,0,1):  -68 y(0,4,1) - u(4,0,1) <= -0
MCLin8(4,0,1):  87 y(0,4,1) + v(4,0,1) - u(4,0,1) <= 87
MCLin9(4,0,1):  87 y(0,4,1) - u(4,0,1) >= -0
MCLin10(4,0,1):  -68 y(0,4,1) + v(4,0,1) - u(4,0,1) >= -68
MCLinV(0,4,2):  x(0,2) - v(0,4,2) - x(4,2) = 2
MCLinV(4,0,2):  - x(0,2) - v(4,0,2) + x(4,2) = 5
MCLin1(0,4,2):  u(0,4,2) >= -0
MCLin2(0,4,2):  -82 y(4,0,2) - u(0,4,2) <= -0
MCLin3(0,4,2):  73 y(4,0,2) + v(0,4,2) - u(0,4,2) <= 73
MCLin4(0,4,2):  73 y(4,0,2) - u(0,4,2) >= -0
MCLin5(0,4,2):  -82 y(4,0,2) + v(0,4,2) - u(0,4,2) >= -82
MCLin6(4,0,2):  u(4,0,2) >= -0
MCLin7(4,0,2):  -80 y(0,4,2) - u(4,0,2) <= -0
MCLin8(4,0,2):  75 y(0,4,2) + v(4,0,2) - u(4,0,2) <= 75
MCLin9(4,0,2):  75 y(0,4,2) - u(4,0,2) >= -0
MCLin10(4,0,2):  -80 y(0,4,2) + v(4,0,2) - u(4,0,2) >= -80
MCLinV(0,4,3):  x(0,3) - v(0,4,3) - x(4,3) = 4
MCLinV(4,0,3):  - x(0,3) - v(4,0,3) + x(4,3) = 7
MCLin1(0,4,3):  u(0,4,3) >= -0
MCLin2(0,4,3):  -74 y(4,0,3) - u(0,4,3) <= -0
MCLin3(0,4,3):  81 y(4,0,3) + v(0,4,3) - u(0,4,3) <= 81
MCLin4(0,4,3):  81 y(4,0,3) - u(0,4,3) >= -0
MCLin5(0,4,3):  -74 y(4,0,3) + v(0,4,3) - u(0,4,3) >= -74
MCLin6(4,0,3):  u(4,0,3) >= -0
MCLin7(4,0,3):  -92 y(0,4,3) - u(4,0,3) <= -0
MCLin8(4,0,3):  63 y(0,4,3) + v(4,0,3) - u(4,0,3) <= 63
MCLin9(4,0,3):  63 y(0,4,3) - u(4,0,3) >= -0
MCLin10(4,0,3):  -92 y(0,4,3) + v(4,0,3) - u(4,0,3) >= -92
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 3
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 4
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -89 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  68 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 68
MCLin4(1,2,0):  68 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -89 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -89
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -75 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  82 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 82
MCLin9(2,1,0):  82 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -75 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -75
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 8
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 2
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -72 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  85 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 85
MCLin4(1,2,1):  85 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -72 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -72
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -95 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  62 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 62
MCLin9(2,1,1):  62 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -95 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -95
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 10
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 2
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -90 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  67 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 67
MCLin4(1,2,2):  67 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -90 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -90
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -79 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  78 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 78
MCLin9(2,1,2):  78 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -79 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -79
MCLinV(1,2,3):  x(1,3) - v(1,2,3) - x(2,3) = 2
MCLinV(2,1,3):  - x(1,3) - v(2,1,3) + x(2,3) = 2
MCLin1(1,2,3):  u(1,2,3) >= -0
MCLin2(1,2,3):  -93 y(2,1,3) - u(1,2,3) <= -0
MCLin3(1,2,3):  64 y(2,1,3) + v(1,2,3) - u(1,2,3) <= 64
MCLin4(1,2,3):  64 y(2,1,3) - u(1,2,3) >= -0
MCLin5(1,2,3):  -93 y(2,1,3) + v(1,2,3) - u(1,2,3) >= -93
MCLin6(2,1,3):  u(2,1,3) >= -0
MCLin7(2,1,3):  -68 y(1,2,3) - u(2,1,3) <= -0
MCLin8(2,1,3):  89 y(1,2,3) + v(2,1,3) - u(2,1,3) <= 89
MCLin9(2,1,3):  89 y(1,2,3) - u(2,1,3) >= -0
MCLin10(2,1,3):  -68 y(1,2,3) + v(2,1,3) - u(2,1,3) >= -68
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 5
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 4
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -91 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  62 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 62
MCLin4(1,3,0):  62 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -91 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -91
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -71 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  82 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 82
MCLin9(3,1,0):  82 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -71 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -71
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 6
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 2
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -75 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  78 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 78
MCLin4(1,3,1):  78 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -75 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -75
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -86 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  67 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 67
MCLin9(3,1,1):  67 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -86 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -86
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 9
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 2
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -77 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  76 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 76
MCLin4(1,3,2):  76 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -77 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -77
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -87 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  66 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 66
MCLin9(3,1,2):  66 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -87 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -87
MCLinV(1,3,3):  x(1,3) - v(1,3,3) - x(3,3) = 7
MCLinV(3,1,3):  - x(1,3) - v(3,1,3) + x(3,3) = 2
MCLin1(1,3,3):  u(1,3,3) >= -0
MCLin2(1,3,3):  -88 y(3,1,3) - u(1,3,3) <= -0
MCLin3(1,3,3):  65 y(3,1,3) + v(1,3,3) - u(1,3,3) <= 65
MCLin4(1,3,3):  65 y(3,1,3) - u(1,3,3) >= -0
MCLin5(1,3,3):  -88 y(3,1,3) + v(1,3,3) - u(1,3,3) >= -88
MCLin6(3,1,3):  u(3,1,3) >= -0
MCLin7(3,1,3):  -74 y(1,3,3) - u(3,1,3) <= -0
MCLin8(3,1,3):  79 y(1,3,3) + v(3,1,3) - u(3,1,3) <= 79
MCLin9(3,1,3):  79 y(1,3,3) - u(3,1,3) >= -0
MCLin10(3,1,3):  -74 y(1,3,3) + v(3,1,3) - u(3,1,3) >= -74
MCLinV(1,4,0):  x(1,0) - v(1,4,0) - x(4,0) = 1
MCLinV(4,1,0):  - x(1,0) - v(4,1,0) + x(4,0) = 4
MCLin1(1,4,0):  u(1,4,0) >= -0
MCLin2(1,4,0):  -84 y(4,1,0) - u(1,4,0) <= -0
MCLin3(1,4,0):  86 y(4,1,0) + v(1,4,0) - u(1,4,0) <= 86
MCLin4(1,4,0):  86 y(4,1,0) - u(1,4,0) >= -0
MCLin5(1,4,0):  -84 y(4,1,0) + v(1,4,0) - u(1,4,0) >= -84
MCLin6(4,1,0):  u(4,1,0) >= -0
MCLin7(4,1,0):  -91 y(1,4,0) - u(4,1,0) <= -0
MCLin8(4,1,0):  79 y(1,4,0) + v(4,1,0) - u(4,1,0) <= 79
MCLin9(4,1,0):  79 y(1,4,0) - u(4,1,0) >= -0
MCLin10(4,1,0):  -91 y(1,4,0) + v(4,1,0) - u(4,1,0) >= -91
MCLinV(1,4,1):  x(1,1) - v(1,4,1) - x(4,1) = 3
MCLinV(4,1,1):  - x(1,1) - v(4,1,1) + x(4,1) = 2
MCLin1(1,4,1):  u(1,4,1) >= -0
MCLin2(1,4,1):  -87 y(4,1,1) - u(1,4,1) <= -0
MCLin3(1,4,1):  83 y(4,1,1) + v(1,4,1) - u(1,4,1) <= 83
MCLin4(1,4,1):  83 y(4,1,1) - u(1,4,1) >= -0
MCLin5(1,4,1):  -87 y(4,1,1) + v(1,4,1) - u(1,4,1) >= -87
MCLin6(4,1,1):  u(4,1,1) >= -0
MCLin7(4,1,1):  -88 y(1,4,1) - u(4,1,1) <= -0
MCLin8(4,1,1):  82 y(1,4,1) + v(4,1,1) - u(4,1,1) <= 82
MCLin9(4,1,1):  82 y(1,4,1) - u(4,1,1) >= -0
MCLin10(4,1,1):  -88 y(1,4,1) + v(4,1,1) - u(4,1,1) >= -88
MCLinV(1,4,2):  x(1,2) - v(1,4,2) - x(4,2) = 2
MCLinV(4,1,2):  - x(1,2) - v(4,1,2) + x(4,2) = 2
MCLin1(1,4,2):  u(1,4,2) >= -0
MCLin2(1,4,2):  -87 y(4,1,2) - u(1,4,2) <= -0
MCLin3(1,4,2):  83 y(4,1,2) + v(1,4,2) - u(1,4,2) <= 83
MCLin4(1,4,2):  83 y(4,1,2) - u(1,4,2) >= -0
MCLin5(1,4,2):  -87 y(4,1,2) + v(1,4,2) - u(1,4,2) >= -87
MCLin6(4,1,2):  u(4,1,2) >= -0
MCLin7(4,1,2):  -87 y(1,4,2) - u(4,1,2) <= -0
MCLin8(4,1,2):  83 y(1,4,2) + v(4,1,2) - u(4,1,2) <= 83
MCLin9(4,1,2):  83 y(1,4,2) - u(4,1,2) >= -0
MCLin10(4,1,2):  -87 y(1,4,2) + v(4,1,2) - u(4,1,2) >= -87
MCLinV(1,4,3):  x(1,3) - v(1,4,3) - x(4,3) = 4
MCLinV(4,1,3):  - x(1,3) - v(4,1,3) + x(4,3) = 2
MCLin1(1,4,3):  u(1,4,3) >= -0
MCLin2(1,4,3):  -90 y(4,1,3) - u(1,4,3) <= -0
MCLin3(1,4,3):  80 y(4,1,3) + v(1,4,3) - u(1,4,3) <= 80
MCLin4(1,4,3):  80 y(4,1,3) - u(1,4,3) >= -0
MCLin5(1,4,3):  -90 y(4,1,3) + v(1,4,3) - u(1,4,3) >= -90
MCLin6(4,1,3):  u(4,1,3) >= -0
MCLin7(4,1,3):  -86 y(1,4,3) - u(4,1,3) <= -0
MCLin8(4,1,3):  84 y(1,4,3) + v(4,1,3) - u(4,1,3) <= 84
MCLin9(4,1,3):  84 y(1,4,3) - u(4,1,3) >= -0
MCLin10(4,1,3):  -86 y(1,4,3) + v(4,1,3) - u(4,1,3) >= -86
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 5
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 3
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -77 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  63 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 63
MCLin4(2,3,0):  63 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -77 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -77
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -71 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  69 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 69
MCLin9(3,2,0):  69 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -71 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -71
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 6
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 8
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -83 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  57 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 57
MCLin4(2,3,1):  57 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -83 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -83
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -71 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  69 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 69
MCLin9(3,2,1):  69 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -71 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -71
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 9
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 10
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -69 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  71 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 71
MCLin4(2,3,2):  71 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -69 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -69
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -90 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  50 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 50
MCLin9(3,2,2):  50 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -90 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -90
MCLinV(2,3,3):  x(2,3) - v(2,3,3) - x(3,3) = 7
MCLinV(3,2,3):  - x(2,3) - v(3,2,3) + x(3,3) = 2
MCLin1(2,3,3):  u(2,3,3) >= -0
MCLin2(2,3,3):  -69 y(3,2,3) - u(2,3,3) <= -0
MCLin3(2,3,3):  71 y(3,2,3) + v(2,3,3) - u(2,3,3) <= 71
MCLin4(2,3,3):  71 y(3,2,3) - u(2,3,3) >= -0
MCLin5(2,3,3):  -69 y(3,2,3) + v(2,3,3) - u(2,3,3) >= -69
MCLin6(3,2,3):  u(3,2,3) >= -0
MCLin7(3,2,3):  -80 y(2,3,3) - u(3,2,3) <= -0
MCLin8(3,2,3):  60 y(2,3,3) + v(3,2,3) - u(3,2,3) <= 60
MCLin9(3,2,3):  60 y(2,3,3) - u(3,2,3) >= -0
MCLin10(3,2,3):  -80 y(2,3,3) + v(3,2,3) - u(3,2,3) >= -80
MCLinV(2,4,0):  x(2,0) - v(2,4,0) - x(4,0) = 1
MCLinV(4,2,0):  - x(2,0) - v(4,2,0) + x(4,0) = 3
MCLin1(2,4,0):  u(2,4,0) >= -0
MCLin2(2,4,0):  -70 y(4,2,0) - u(2,4,0) <= -0
MCLin3(2,4,0):  87 y(4,2,0) + v(2,4,0) - u(2,4,0) <= 87
MCLin4(2,4,0):  87 y(4,2,0) - u(2,4,0) >= -0
MCLin5(2,4,0):  -70 y(4,2,0) + v(2,4,0) - u(2,4,0) >= -70
MCLin6(4,2,0):  u(4,2,0) >= -0
MCLin7(4,2,0):  -91 y(2,4,0) - u(4,2,0) <= -0
MCLin8(4,2,0):  66 y(2,4,0) + v(4,2,0) - u(4,2,0) <= 66
MCLin9(4,2,0):  66 y(2,4,0) - u(4,2,0) >= -0
MCLin10(4,2,0):  -91 y(2,4,0) + v(4,2,0) - u(4,2,0) >= -91
MCLinV(2,4,1):  x(2,1) - v(2,4,1) - x(4,1) = 3
MCLinV(4,2,1):  - x(2,1) - v(4,2,1) + x(4,1) = 8
MCLin1(2,4,1):  u(2,4,1) >= -0
MCLin2(2,4,1):  -95 y(4,2,1) - u(2,4,1) <= -0
MCLin3(2,4,1):  62 y(4,2,1) + v(2,4,1) - u(2,4,1) <= 62
MCLin4(2,4,1):  62 y(4,2,1) - u(2,4,1) >= -0
MCLin5(2,4,1):  -95 y(4,2,1) + v(2,4,1) - u(2,4,1) >= -95
MCLin6(4,2,1):  u(4,2,1) >= -0
MCLin7(4,2,1):  -73 y(2,4,1) - u(4,2,1) <= -0
MCLin8(4,2,1):  84 y(2,4,1) + v(4,2,1) - u(4,2,1) <= 84
MCLin9(4,2,1):  84 y(2,4,1) - u(4,2,1) >= -0
MCLin10(4,2,1):  -73 y(2,4,1) + v(4,2,1) - u(4,2,1) >= -73
MCLinV(2,4,2):  x(2,2) - v(2,4,2) - x(4,2) = 2
MCLinV(4,2,2):  - x(2,2) - v(4,2,2) + x(4,2) = 10
MCLin1(2,4,2):  u(2,4,2) >= -0
MCLin2(2,4,2):  -79 y(4,2,2) - u(2,4,2) <= -0
MCLin3(2,4,2):  78 y(4,2,2) + v(2,4,2) - u(2,4,2) <= 78
MCLin4(2,4,2):  78 y(4,2,2) - u(2,4,2) >= -0
MCLin5(2,4,2):  -79 y(4,2,2) + v(2,4,2) - u(2,4,2) >= -79
MCLin6(4,2,2):  u(4,2,2) >= -0
MCLin7(4,2,2):  -90 y(2,4,2) - u(4,2,2) <= -0
MCLin8(4,2,2):  67 y(2,4,2) + v(4,2,2) - u(4,2,2) <= 67
MCLin9(4,2,2):  67 y(2,4,2) - u(4,2,2) >= -0
MCLin10(4,2,2):  -90 y(2,4,2) + v(4,2,2) - u(4,2,2) >= -90
MCLinV(2,4,3):  x(2,3) - v(2,4,3) - x(4,3) = 4
MCLinV(4,2,3):  - x(2,3) - v(4,2,3) + x(4,3) = 2
MCLin1(2,4,3):  u(2,4,3) >= -0
MCLin2(2,4,3):  -71 y(4,2,3) - u(2,4,3) <= -0
MCLin3(2,4,3):  86 y(4,2,3) + v(2,4,3) - u(2,4,3) <= 86
MCLin4(2,4,3):  86 y(4,2,3) - u(2,4,3) >= -0
MCLin5(2,4,3):  -71 y(4,2,3) + v(2,4,3) - u(2,4,3) >= -71
MCLin6(4,2,3):  u(4,2,3) >= -0
MCLin7(4,2,3):  -92 y(2,4,3) - u(4,2,3) <= -0
MCLin8(4,2,3):  65 y(2,4,3) + v(4,2,3) - u(4,2,3) <= 65
MCLin9(4,2,3):  65 y(2,4,3) - u(4,2,3) >= -0
MCLin10(4,2,3):  -92 y(2,4,3) + v(4,2,3) - u(4,2,3) >= -92
MCLinV(3,4,0):  x(3,0) - v(3,4,0) - x(4,0) = 1
MCLinV(4,3,0):  - x(3,0) - v(4,3,0) + x(4,0) = 5
MCLin1(3,4,0):  u(3,4,0) >= -0
MCLin2(3,4,0):  -66 y(4,3,0) - u(3,4,0) <= -0
MCLin3(3,4,0):  87 y(4,3,0) + v(3,4,0) - u(3,4,0) <= 87
MCLin4(3,4,0):  87 y(4,3,0) - u(3,4,0) >= -0
MCLin5(3,4,0):  -66 y(4,3,0) + v(3,4,0) - u(3,4,0) >= -66
MCLin6(4,3,0):  u(4,3,0) >= -0
MCLin7(4,3,0):  -93 y(3,4,0) - u(4,3,0) <= -0
MCLin8(4,3,0):  60 y(3,4,0) + v(4,3,0) - u(4,3,0) <= 60
MCLin9(4,3,0):  60 y(3,4,0) - u(4,3,0) >= -0
MCLin10(4,3,0):  -93 y(3,4,0) + v(4,3,0) - u(4,3,0) >= -93
MCLinV(3,4,1):  x(3,1) - v(3,4,1) - x(4,1) = 3
MCLinV(4,3,1):  - x(3,1) - v(4,3,1) + x(4,1) = 6
MCLin1(3,4,1):  u(3,4,1) >= -0
MCLin2(3,4,1):  -86 y(4,3,1) - u(3,4,1) <= -0
MCLin3(3,4,1):  67 y(4,3,1) + v(3,4,1) - u(3,4,1) <= 67
MCLin4(3,4,1):  67 y(4,3,1) - u(3,4,1) >= -0
MCLin5(3,4,1):  -86 y(4,3,1) + v(3,4,1) - u(3,4,1) >= -86
MCLin6(4,3,1):  u(4,3,1) >= -0
MCLin7(4,3,1):  -76 y(3,4,1) - u(4,3,1) <= -0
MCLin8(4,3,1):  77 y(3,4,1) + v(4,3,1) - u(4,3,1) <= 77
MCLin9(4,3,1):  77 y(3,4,1) - u(4,3,1) >= -0
MCLin10(4,3,1):  -76 y(3,4,1) + v(4,3,1) - u(4,3,1) >= -76
MCLinV(3,4,2):  x(3,2) - v(3,4,2) - x(4,2) = 2
MCLinV(4,3,2):  - x(3,2) - v(4,3,2) + x(4,2) = 9
MCLin1(3,4,2):  u(3,4,2) >= -0
MCLin2(3,4,2):  -87 y(4,3,2) - u(3,4,2) <= -0
MCLin3(3,4,2):  66 y(4,3,2) + v(3,4,2) - u(3,4,2) <= 66
MCLin4(3,4,2):  66 y(4,3,2) - u(3,4,2) >= -0
MCLin5(3,4,2):  -87 y(4,3,2) + v(3,4,2) - u(3,4,2) >= -87
MCLin6(4,3,2):  u(4,3,2) >= -0
MCLin7(4,3,2):  -77 y(3,4,2) - u(4,3,2) <= -0
MCLin8(4,3,2):  76 y(3,4,2) + v(4,3,2) - u(4,3,2) <= 76
MCLin9(4,3,2):  76 y(3,4,2) - u(4,3,2) >= -0
MCLin10(4,3,2):  -77 y(3,4,2) + v(4,3,2) - u(4,3,2) >= -77
MCLinV(3,4,3):  x(3,3) - v(3,4,3) - x(4,3) = 4
MCLinV(4,3,3):  - x(3,3) - v(4,3,3) + x(4,3) = 7
MCLin1(3,4,3):  u(3,4,3) >= -0
MCLin2(3,4,3):  -77 y(4,3,3) - u(3,4,3) <= -0
MCLin3(3,4,3):  76 y(4,3,3) + v(3,4,3) - u(3,4,3) <= 76
MCLin4(3,4,3):  76 y(4,3,3) - u(3,4,3) >= -0
MCLin5(3,4,3):  -77 y(4,3,3) + v(3,4,3) - u(3,4,3) >= -77
MCLin6(4,3,3):  u(4,3,3) >= -0
MCLin7(4,3,3):  -87 y(3,4,3) - u(4,3,3) <= -0
MCLin8(4,3,3):  66 y(3,4,3) + v(4,3,3) - u(4,3,3) <= 66
MCLin9(4,3,3):  66 y(3,4,3) - u(4,3,3) >= -0
MCLin10(4,3,3):  -87 y(3,4,3) + v(4,3,3) - u(4,3,3) >= -87
makespan(0):  C - x(0,3) >= 7
makespan(1):  C - x(1,1) >= 2
makespan(2):  C - x(2,3) >= 2
makespan(3):  C - x(3,0) >= 5
makespan(4):  C - x(4,1) >= 3
Bounds
 0 <= C <= 95
 10 <= x(0,0) <= 80
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -83 <= v(0,1,0) <= 72
 -84 <= v(1,0,0) <= 71
 -950 <= u(0,1,0)
 -950 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -83 <= v(0,2,0) <= 59
 -70 <= v(2,0,0) <= 72
 -950 <= u(0,2,0)
 -950 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -85 <= v(0,3,0) <= 53
 -66 <= v(3,0,0) <= 72
 -950 <= u(0,3,0)
 -950 <= u(3,0,0)
 0 <= y(0,4,0) <= 1
 0 <= y(4,0,0) <= 1
 -78 <= v(0,4,0) <= 77
 -86 <= v(4,0,0) <= 69
 -950 <= u(0,4,0)
 -950 <= u(4,0,0)
 4 <= x(1,0) <= 89
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -89 <= v(1,2,0) <= 68
 -75 <= v(2,1,0) <= 82
 -950 <= u(1,2,0)
 -950 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -91 <= v(1,3,0) <= 62
 -71 <= v(3,1,0) <= 82
 -950 <= u(1,3,0)
 -950 <= u(3,1,0)
 0 <= y(1,4,0) <= 1
 0 <= y(4,1,0) <= 1
 -84 <= v(1,4,0) <= 86
 -91 <= v(4,1,0) <= 79
 -950 <= u(1,4,0)
 -950 <= u(4,1,0)
 18 <= x(2,0) <= 90
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -77 <= v(2,3,0) <= 63
 -71 <= v(3,2,0) <= 69
 -950 <= u(2,3,0)
 -950 <= u(3,2,0)
 0 <= y(2,4,0) <= 1
 0 <= y(4,2,0) <= 1
 -70 <= v(2,4,0) <= 87
 -91 <= v(4,2,0) <= 66
 -950 <= u(2,4,0)
 -950 <= u(4,2,0)
 22 <= x(3,0) <= 90
 0 <= y(3,4,0) <= 1
 0 <= y(4,3,0) <= 1
 -66 <= v(3,4,0) <= 87
 -93 <= v(4,3,0) <= 60
 -950 <= u(3,4,0)
 -950 <= u(4,3,0)
 2 <= x(4,0) <= 87
 0 <= x(0,1) <= 70
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -95 <= v(0,1,1) <= 60
 -67 <= v(1,0,1) <= 88
 -950 <= u(0,1,1)
 -950 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -80 <= v(0,2,1) <= 62
 -75 <= v(2,0,1) <= 67
 -950 <= u(0,2,1)
 -950 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -83 <= v(0,3,1) <= 55
 -66 <= v(3,0,1) <= 72
 -950 <= u(0,3,1)
 -950 <= u(3,0,1)
 0 <= y(0,4,1) <= 1
 0 <= y(4,0,1) <= 1
 -95 <= v(0,4,1) <= 60
 -68 <= v(4,0,1) <= 87
 -950 <= u(0,4,1)
 -950 <= u(4,0,1)
 8 <= x(1,1) <= 93
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -72 <= v(1,2,1) <= 85
 -95 <= v(2,1,1) <= 62
 -950 <= u(1,2,1)
 -950 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -75 <= v(1,3,1) <= 78
 -86 <= v(3,1,1) <= 67
 -950 <= u(1,3,1)
 -950 <= u(3,1,1)
 0 <= y(1,4,1) <= 1
 0 <= y(4,1,1) <= 1
 -87 <= v(1,4,1) <= 83
 -88 <= v(4,1,1) <= 82
 -950 <= u(1,4,1)
 -950 <= u(4,1,1)
 0 <= x(2,1) <= 72
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -83 <= v(2,3,1) <= 57
 -71 <= v(3,2,1) <= 69
 -950 <= u(2,3,1)
 -950 <= u(3,2,1)
 0 <= y(2,4,1) <= 1
 0 <= y(4,2,1) <= 1
 -95 <= v(2,4,1) <= 62
 -73 <= v(4,2,1) <= 84
 -950 <= u(2,4,1)
 -950 <= u(4,2,1)
 9 <= x(3,1) <= 77
 0 <= y(3,4,1) <= 1
 0 <= y(4,3,1) <= 1
 -86 <= v(3,4,1) <= 67
 -76 <= v(4,3,1) <= 77
 -950 <= u(3,4,1)
 -950 <= u(4,3,1)
 7 <= x(4,1) <= 92
 5 <= x(0,2) <= 75
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -82 <= v(0,1,2) <= 73
 -80 <= v(1,0,2) <= 75
 -950 <= u(0,1,2)
 -950 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -85 <= v(0,2,2) <= 57
 -72 <= v(2,0,2) <= 70
 -950 <= u(0,2,2)
 -950 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -72 <= v(0,3,2) <= 66
 -80 <= v(3,0,2) <= 58
 -950 <= u(0,3,2)
 -950 <= u(3,0,2)
 0 <= y(0,4,2) <= 1
 0 <= y(4,0,2) <= 1
 -82 <= v(0,4,2) <= 73
 -80 <= v(4,0,2) <= 75
 -950 <= u(0,4,2)
 -950 <= u(4,0,2)
 0 <= x(1,2) <= 85
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -90 <= v(1,2,2) <= 67
 -79 <= v(2,1,2) <= 78
 -950 <= u(1,2,2)
 -950 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -77 <= v(1,3,2) <= 76
 -87 <= v(3,1,2) <= 66
 -950 <= u(1,3,2)
 -950 <= u(3,1,2)
 0 <= y(1,4,2) <= 1
 0 <= y(4,1,2) <= 1
 -87 <= v(1,4,2) <= 83
 -87 <= v(4,1,2) <= 83
 -950 <= u(1,4,2)
 -950 <= u(4,1,2)
 8 <= x(2,2) <= 80
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -69 <= v(2,3,2) <= 71
 -90 <= v(3,2,2) <= 50
 -950 <= u(2,3,2)
 -950 <= u(3,2,2)
 0 <= y(2,4,2) <= 1
 0 <= y(4,2,2) <= 1
 -79 <= v(2,4,2) <= 78
 -90 <= v(4,2,2) <= 67
 -950 <= u(2,4,2)
 -950 <= u(4,2,2)
 0 <= x(3,2) <= 68
 0 <= y(3,4,2) <= 1
 0 <= y(4,3,2) <= 1
 -87 <= v(3,4,2) <= 66
 -77 <= v(4,3,2) <= 76
 -950 <= u(3,4,2)
 -950 <= u(4,3,2)
 0 <= x(4,2) <= 85
 18 <= x(0,3) <= 88
 0 <= y(0,1,3) <= 1
 0 <= y(1,0,3) <= 1
 -71 <= v(0,1,3) <= 84
 -93 <= v(1,0,3) <= 62
 -950 <= u(0,1,3)
 -950 <= u(1,0,3)
 0 <= y(0,2,3) <= 1
 0 <= y(2,0,3) <= 1
 -77 <= v(0,2,3) <= 65
 -74 <= v(2,0,3) <= 68
 -950 <= u(0,2,3)
 -950 <= u(2,0,3)
 0 <= y(0,3,3) <= 1
 0 <= y(3,0,3) <= 1
 -72 <= v(0,3,3) <= 66
 -80 <= v(3,0,3) <= 58
 -950 <= u(0,3,3)
 -950 <= u(3,0,3)
 0 <= y(0,4,3) <= 1
 0 <= y(4,0,3) <= 1
 -74 <= v(0,4,3) <= 81
 -92 <= v(4,0,3) <= 63
 -950 <= u(0,4,3)
 -950 <= u(4,0,3)
 2 <= x(1,3) <= 87
 0 <= y(1,2,3) <= 1
 0 <= y(2,1,3) <= 1
 -93 <= v(1,2,3) <= 64
 -68 <= v(2,1,3) <= 89
 -950 <= u(1,2,3)
 -950 <= u(2,1,3)
 0 <= y(1,3,3) <= 1
 0 <= y(3,1,3) <= 1
 -88 <= v(1,3,3) <= 65
 -74 <= v(3,1,3) <= 79
 -950 <= u(1,3,3)
 -950 <= u(3,1,3)
 0 <= y(1,4,3) <= 1
 0 <= y(4,1,3) <= 1
 -90 <= v(1,4,3) <= 80
 -86 <= v(4,1,3) <= 84
 -950 <= u(1,4,3)
 -950 <= u(4,1,3)
 21 <= x(2,3) <= 93
 0 <= y(2,3,3) <= 1
 0 <= y(3,2,3) <= 1
 -69 <= v(2,3,3) <= 71
 -80 <= v(3,2,3) <= 60
 -950 <= u(2,3,3)
 -950 <= u(3,2,3)
 0 <= y(2,4,3) <= 1
 0 <= y(4,2,3) <= 1
 -71 <= v(2,4,3) <= 86
 -92 <= v(4,2,3) <= 65
 -950 <= u(2,4,3)
 -950 <= u(4,2,3)
 15 <= x(3,3) <= 83
 0 <= y(3,4,3) <= 1
 0 <= y(4,3,3) <= 1
 -77 <= v(3,4,3) <= 76
 -87 <= v(4,3,3) <= 66
 -950 <= u(3,4,3)
 -950 <= u(4,3,3)
 3 <= x(4,3) <= 88
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
y(4,3,2) v(3,4,2) v(4,3,2) u(3,4,2) u(4,3,2) x(4,2) x(0,3) y(0,1,3) y(1,0,3) v(0,1,3) 
v(1,0,3) u(0,1,3) u(1,0,3) y(0,2,3) y(2,0,3) v(0,2,3) v(2,0,3) u(0,2,3) u(2,0,3) y(0,3,3) 
y(3,0,3) v(0,3,3) v(3,0,3) u(0,3,3) u(3,0,3) y(0,4,3) y(4,0,3) v(0,4,3) v(4,0,3) u(0,4,3) 
u(4,0,3) x(1,3) y(1,2,3) y(2,1,3) v(1,2,3) v(2,1,3) u(1,2,3) u(2,1,3) y(1,3,3) y(3,1,3) 
v(1,3,3) v(3,1,3) u(1,3,3) u(3,1,3) y(1,4,3) y(4,1,3) v(1,4,3) v(4,1,3) u(1,4,3) u(4,1,3) 
x(2,3) y(2,3,3) y(3,2,3) v(2,3,3) v(3,2,3) u(2,3,3) u(3,2,3) y(2,4,3) y(4,2,3) v(2,4,3) 
v(4,2,3) u(2,4,3) u(4,2,3) x(3,3) y(3,4,3) y(4,3,3) v(3,4,3) v(4,3,3) u(3,4,3) u(4,3,3) 
x(4,3) 
End
