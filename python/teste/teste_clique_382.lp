\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,2) >= 5
ord(0,2):  x(0,0) - x(0,2) >= 7
ord(0,3):  - x(0,0) + x(0,3) >= 8
ord(1,1):  - x(1,0) + x(1,1) >= 4
ord(1,2):  - x(1,1) + x(1,3) >= 3
ord(1,3):  x(1,2) - x(1,3) >= 10
ord(2,1):  x(2,1) - x(2,2) >= 9
ord(2,2):  - x(2,1) + x(2,3) >= 2
ord(2,3):  x(2,0) - x(2,3) >= 9
ord(3,1):  - x(3,1) + x(3,3) >= 1
ord(3,2):  x(3,0) - x(3,3) >= 2
ord(3,3):  - x(3,0) + x(3,2) >= 3
ord(4,1):  x(4,0) - x(4,3) >= 1
ord(4,2):  - x(4,0) + x(4,2) >= 8
ord(4,3):  x(4,1) - x(4,2) >= 1
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
MCLin2(0,1,0):  -62 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  76 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 76
MCLin4(0,1,0):  76 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -62 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -62
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -88 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  50 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 50
MCLin9(1,0,0):  50 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -88 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -88
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 3
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 5
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -77 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  61 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 61
MCLin4(0,1,1):  61 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -77 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -77
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -69 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  69 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 69
MCLin9(1,0,1):  69 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -69 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -69
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 2
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 7
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -84 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  54 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 54
MCLin4(0,1,2):  54 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -84 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -84
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -63 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  75 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 75
MCLin9(1,0,2):  75 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -63 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -63
MCLinV(0,1,3):  x(0,3) - v(0,1,3) - x(1,3) = 10
MCLinV(1,0,3):  - x(0,3) - v(1,0,3) + x(1,3) = 1
MCLin1(0,1,3):  u(0,1,3) >= -0
MCLin2(0,1,3):  -67 y(1,0,3) - u(0,1,3) <= -0
MCLin3(0,1,3):  71 y(1,0,3) + v(0,1,3) - u(0,1,3) <= 71
MCLin4(0,1,3):  71 y(1,0,3) - u(0,1,3) >= -0
MCLin5(0,1,3):  -67 y(1,0,3) + v(0,1,3) - u(0,1,3) >= -67
MCLin6(1,0,3):  u(1,0,3) >= -0
MCLin7(1,0,3):  -82 y(0,1,3) - u(1,0,3) <= -0
MCLin8(1,0,3):  56 y(0,1,3) + v(1,0,3) - u(1,0,3) <= 56
MCLin9(1,0,3):  56 y(0,1,3) - u(1,0,3) >= -0
MCLin10(1,0,3):  -82 y(0,1,3) + v(1,0,3) - u(1,0,3) >= -82
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 3
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 8
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -77 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  57 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 57
MCLin4(0,2,0):  57 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -77 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -77
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -68 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  66 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 66
MCLin9(2,0,0):  66 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -68 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -68
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 5
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -77 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  57 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 57
MCLin4(0,2,1):  57 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -77 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -77
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -64 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  70 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 70
MCLin9(2,0,1):  70 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -64 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -64
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 9
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 7
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -70 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  64 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 64
MCLin4(0,2,2):  64 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -70 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -70
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -80 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  54 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 54
MCLin9(2,0,2):  54 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -80 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -80
MCLinV(0,2,3):  x(0,3) - v(0,2,3) - x(2,3) = 9
MCLinV(2,0,3):  - x(0,3) - v(2,0,3) + x(2,3) = 1
MCLin1(0,2,3):  u(0,2,3) >= -0
MCLin2(0,2,3):  -66 y(2,0,3) - u(0,2,3) <= -0
MCLin3(0,2,3):  68 y(2,0,3) + v(0,2,3) - u(0,2,3) <= 68
MCLin4(0,2,3):  68 y(2,0,3) - u(0,2,3) >= -0
MCLin5(0,2,3):  -66 y(2,0,3) + v(0,2,3) - u(0,2,3) >= -66
MCLin6(2,0,3):  u(2,0,3) >= -0
MCLin7(2,0,3):  -78 y(0,2,3) - u(2,0,3) <= -0
MCLin8(2,0,3):  56 y(0,2,3) + v(2,0,3) - u(2,0,3) <= 56
MCLin9(2,0,3):  56 y(0,2,3) - u(2,0,3) >= -0
MCLin10(2,0,3):  -78 y(0,2,3) + v(2,0,3) - u(2,0,3) >= -78
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 3
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 8
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -76 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  74 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 74
MCLin4(0,3,0):  74 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -76 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -76
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -85 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  65 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 65
MCLin9(3,0,0):  65 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -85 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -85
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 1
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 5
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -83 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  67 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 67
MCLin4(0,3,1):  67 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -83 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -83
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -73 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  77 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 77
MCLin9(3,0,1):  77 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -73 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -73
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 1
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 7
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -84 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  66 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 66
MCLin4(0,3,2):  66 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -84 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -84
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -74 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  76 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 76
MCLin9(3,0,2):  76 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -74 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -74
MCLinV(0,3,3):  x(0,3) - v(0,3,3) - x(3,3) = 2
MCLinV(3,0,3):  - x(0,3) - v(3,0,3) + x(3,3) = 1
MCLin1(0,3,3):  u(0,3,3) >= -0
MCLin2(0,3,3):  -65 y(3,0,3) - u(0,3,3) <= -0
MCLin3(0,3,3):  85 y(3,0,3) + v(0,3,3) - u(0,3,3) <= 85
MCLin4(0,3,3):  85 y(3,0,3) - u(0,3,3) >= -0
MCLin5(0,3,3):  -65 y(3,0,3) + v(0,3,3) - u(0,3,3) >= -65
MCLin6(3,0,3):  u(3,0,3) >= -0
MCLin7(3,0,3):  -88 y(0,3,3) - u(3,0,3) <= -0
MCLin8(3,0,3):  62 y(0,3,3) + v(3,0,3) - u(3,0,3) <= 62
MCLin9(3,0,3):  62 y(0,3,3) - u(3,0,3) >= -0
MCLin10(3,0,3):  -88 y(0,3,3) + v(3,0,3) - u(3,0,3) >= -88
MCLinV(0,4,0):  x(0,0) - v(0,4,0) - x(4,0) = 8
MCLinV(4,0,0):  - x(0,0) - v(4,0,0) + x(4,0) = 8
MCLin1(0,4,0):  u(0,4,0) >= -0
MCLin2(0,4,0):  -67 y(4,0,0) - u(0,4,0) <= -0
MCLin3(0,4,0):  71 y(4,0,0) + v(0,4,0) - u(0,4,0) <= 71
MCLin4(0,4,0):  71 y(4,0,0) - u(0,4,0) >= -0
MCLin5(0,4,0):  -67 y(4,0,0) + v(0,4,0) - u(0,4,0) >= -67
MCLin6(4,0,0):  u(4,0,0) >= -0
MCLin7(4,0,0):  -87 y(0,4,0) - u(4,0,0) <= -0
MCLin8(4,0,0):  51 y(0,4,0) + v(4,0,0) - u(4,0,0) <= 51
MCLin9(4,0,0):  51 y(0,4,0) - u(4,0,0) >= -0
MCLin10(4,0,0):  -87 y(0,4,0) + v(4,0,0) - u(4,0,0) >= -87
MCLinV(0,4,1):  x(0,1) - v(0,4,1) - x(4,1) = 9
MCLinV(4,0,1):  - x(0,1) - v(4,0,1) + x(4,1) = 5
MCLin1(0,4,1):  u(0,4,1) >= -0
MCLin2(0,4,1):  -89 y(4,0,1) - u(0,4,1) <= -0
MCLin3(0,4,1):  49 y(4,0,1) + v(0,4,1) - u(0,4,1) <= 49
MCLin4(0,4,1):  49 y(4,0,1) - u(0,4,1) >= -0
MCLin5(0,4,1):  -89 y(4,0,1) + v(0,4,1) - u(0,4,1) >= -89
MCLin6(4,0,1):  u(4,0,1) >= -0
MCLin7(4,0,1):  -63 y(0,4,1) - u(4,0,1) <= -0
MCLin8(4,0,1):  75 y(0,4,1) + v(4,0,1) - u(4,0,1) <= 75
MCLin9(4,0,1):  75 y(0,4,1) - u(4,0,1) >= -0
MCLin10(4,0,1):  -63 y(0,4,1) + v(4,0,1) - u(4,0,1) >= -63
MCLinV(0,4,2):  x(0,2) - v(0,4,2) - x(4,2) = 1
MCLinV(4,0,2):  - x(0,2) - v(4,0,2) + x(4,2) = 7
MCLin1(0,4,2):  u(0,4,2) >= -0
MCLin2(0,4,2):  -75 y(4,0,2) - u(0,4,2) <= -0
MCLin3(0,4,2):  63 y(4,0,2) + v(0,4,2) - u(0,4,2) <= 63
MCLin4(0,4,2):  63 y(4,0,2) - u(0,4,2) >= -0
MCLin5(0,4,2):  -75 y(4,0,2) + v(0,4,2) - u(0,4,2) >= -75
MCLin6(4,0,2):  u(4,0,2) >= -0
MCLin7(4,0,2):  -71 y(0,4,2) - u(4,0,2) <= -0
MCLin8(4,0,2):  67 y(0,4,2) + v(4,0,2) - u(4,0,2) <= 67
MCLin9(4,0,2):  67 y(0,4,2) - u(4,0,2) >= -0
MCLin10(4,0,2):  -71 y(0,4,2) + v(4,0,2) - u(4,0,2) >= -71
MCLinV(0,4,3):  x(0,3) - v(0,4,3) - x(4,3) = 1
MCLinV(4,0,3):  - x(0,3) - v(4,0,3) + x(4,3) = 1
MCLin1(0,4,3):  u(0,4,3) >= -0
MCLin2(0,4,3):  -51 y(4,0,3) - u(0,4,3) <= -0
MCLin3(0,4,3):  87 y(4,0,3) + v(0,4,3) - u(0,4,3) <= 87
MCLin4(0,4,3):  87 y(4,0,3) - u(0,4,3) >= -0
MCLin5(0,4,3):  -51 y(4,0,3) + v(0,4,3) - u(0,4,3) >= -51
MCLin6(4,0,3):  u(4,0,3) >= -0
MCLin7(4,0,3):  -89 y(0,4,3) - u(4,0,3) <= -0
MCLin8(4,0,3):  49 y(0,4,3) + v(4,0,3) - u(4,0,3) <= 49
MCLin9(4,0,3):  49 y(0,4,3) - u(4,0,3) >= -0
MCLin10(4,0,3):  -89 y(0,4,3) + v(4,0,3) - u(4,0,3) >= -89
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 3
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 4
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -89 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  47 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 47
MCLin4(1,2,0):  47 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -89 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -89
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -54 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  82 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 82
MCLin9(2,1,0):  82 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -54 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -54
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 3
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -73 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  63 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 63
MCLin4(1,2,1):  63 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -73 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -73
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -68 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  68 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 68
MCLin9(2,1,1):  68 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -68 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -68
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 9
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 2
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -58 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  78 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 78
MCLin4(1,2,2):  78 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -58 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -58
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -89 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  47 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 47
MCLin9(2,1,2):  47 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -89 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -89
MCLinV(1,2,3):  x(1,3) - v(1,2,3) - x(2,3) = 9
MCLinV(2,1,3):  - x(1,3) - v(2,1,3) + x(2,3) = 10
MCLin1(1,2,3):  u(1,2,3) >= -0
MCLin2(1,2,3):  -79 y(2,1,3) - u(1,2,3) <= -0
MCLin3(1,2,3):  57 y(2,1,3) + v(1,2,3) - u(1,2,3) <= 57
MCLin4(1,2,3):  57 y(2,1,3) - u(1,2,3) >= -0
MCLin5(1,2,3):  -79 y(2,1,3) + v(1,2,3) - u(1,2,3) >= -79
MCLin6(2,1,3):  u(2,1,3) >= -0
MCLin7(2,1,3):  -76 y(1,2,3) - u(2,1,3) <= -0
MCLin8(2,1,3):  60 y(1,2,3) + v(2,1,3) - u(2,1,3) <= 60
MCLin9(2,1,3):  60 y(1,2,3) - u(2,1,3) >= -0
MCLin10(2,1,3):  -76 y(1,2,3) + v(2,1,3) - u(2,1,3) >= -76
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 3
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 4
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -88 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  64 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 64
MCLin4(1,3,0):  64 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -88 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -88
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -71 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  81 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 81
MCLin9(3,1,0):  81 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -71 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -71
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 1
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 3
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -79 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  73 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 73
MCLin4(1,3,1):  73 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -79 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -79
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -77 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  75 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 75
MCLin9(3,1,1):  75 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -77 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -77
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 1
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 2
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -72 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  80 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 80
MCLin4(1,3,2):  80 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -72 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -72
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -83 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  69 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 69
MCLin9(3,1,2):  69 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -83 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -83
MCLinV(1,3,3):  x(1,3) - v(1,3,3) - x(3,3) = 2
MCLinV(3,1,3):  - x(1,3) - v(3,1,3) + x(3,3) = 10
MCLin1(1,3,3):  u(1,3,3) >= -0
MCLin2(1,3,3):  -78 y(3,1,3) - u(1,3,3) <= -0
MCLin3(1,3,3):  74 y(3,1,3) + v(1,3,3) - u(1,3,3) <= 74
MCLin4(1,3,3):  74 y(3,1,3) - u(1,3,3) >= -0
MCLin5(1,3,3):  -78 y(3,1,3) + v(1,3,3) - u(1,3,3) >= -78
MCLin6(3,1,3):  u(3,1,3) >= -0
MCLin7(3,1,3):  -86 y(1,3,3) - u(3,1,3) <= -0
MCLin8(3,1,3):  66 y(1,3,3) + v(3,1,3) - u(3,1,3) <= 66
MCLin9(3,1,3):  66 y(1,3,3) - u(3,1,3) >= -0
MCLin10(3,1,3):  -86 y(1,3,3) + v(3,1,3) - u(3,1,3) >= -86
MCLinV(1,4,0):  x(1,0) - v(1,4,0) - x(4,0) = 8
MCLinV(4,1,0):  - x(1,0) - v(4,1,0) + x(4,0) = 4
MCLin1(1,4,0):  u(1,4,0) >= -0
MCLin2(1,4,0):  -79 y(4,1,0) - u(1,4,0) <= -0
MCLin3(1,4,0):  61 y(4,1,0) + v(1,4,0) - u(1,4,0) <= 61
MCLin4(1,4,0):  61 y(4,1,0) - u(1,4,0) >= -0
MCLin5(1,4,0):  -79 y(4,1,0) + v(1,4,0) - u(1,4,0) >= -79
MCLin6(4,1,0):  u(4,1,0) >= -0
MCLin7(4,1,0):  -73 y(1,4,0) - u(4,1,0) <= -0
MCLin8(4,1,0):  67 y(1,4,0) + v(4,1,0) - u(4,1,0) <= 67
MCLin9(4,1,0):  67 y(1,4,0) - u(4,1,0) >= -0
MCLin10(4,1,0):  -73 y(1,4,0) + v(4,1,0) - u(4,1,0) >= -73
MCLinV(1,4,1):  x(1,1) - v(1,4,1) - x(4,1) = 9
MCLinV(4,1,1):  - x(1,1) - v(4,1,1) + x(4,1) = 3
MCLin1(1,4,1):  u(1,4,1) >= -0
MCLin2(1,4,1):  -85 y(4,1,1) - u(1,4,1) <= -0
MCLin3(1,4,1):  55 y(4,1,1) + v(1,4,1) - u(1,4,1) <= 55
MCLin4(1,4,1):  55 y(4,1,1) - u(1,4,1) >= -0
MCLin5(1,4,1):  -85 y(4,1,1) + v(1,4,1) - u(1,4,1) >= -85
MCLin6(4,1,1):  u(4,1,1) >= -0
MCLin7(4,1,1):  -67 y(1,4,1) - u(4,1,1) <= -0
MCLin8(4,1,1):  73 y(1,4,1) + v(4,1,1) - u(4,1,1) <= 73
MCLin9(4,1,1):  73 y(1,4,1) - u(4,1,1) >= -0
MCLin10(4,1,1):  -67 y(1,4,1) + v(4,1,1) - u(4,1,1) >= -67
MCLinV(1,4,2):  x(1,2) - v(1,4,2) - x(4,2) = 1
MCLinV(4,1,2):  - x(1,2) - v(4,1,2) + x(4,2) = 2
MCLin1(1,4,2):  u(1,4,2) >= -0
MCLin2(1,4,2):  -63 y(4,1,2) - u(1,4,2) <= -0
MCLin3(1,4,2):  77 y(4,1,2) + v(1,4,2) - u(1,4,2) <= 77
MCLin4(1,4,2):  77 y(4,1,2) - u(1,4,2) >= -0
MCLin5(1,4,2):  -63 y(4,1,2) + v(1,4,2) - u(1,4,2) >= -63
MCLin6(4,1,2):  u(4,1,2) >= -0
MCLin7(4,1,2):  -80 y(1,4,2) - u(4,1,2) <= -0
MCLin8(4,1,2):  60 y(1,4,2) + v(4,1,2) - u(4,1,2) <= 60
MCLin9(4,1,2):  60 y(1,4,2) - u(4,1,2) >= -0
MCLin10(4,1,2):  -80 y(1,4,2) + v(4,1,2) - u(4,1,2) >= -80
MCLinV(1,4,3):  x(1,3) - v(1,4,3) - x(4,3) = 1
MCLinV(4,1,3):  - x(1,3) - v(4,1,3) + x(4,3) = 10
MCLin1(1,4,3):  u(1,4,3) >= -0
MCLin2(1,4,3):  -64 y(4,1,3) - u(1,4,3) <= -0
MCLin3(1,4,3):  76 y(4,1,3) + v(1,4,3) - u(1,4,3) <= 76
MCLin4(1,4,3):  76 y(4,1,3) - u(1,4,3) >= -0
MCLin5(1,4,3):  -64 y(4,1,3) + v(1,4,3) - u(1,4,3) >= -64
MCLin6(4,1,3):  u(4,1,3) >= -0
MCLin7(4,1,3):  -87 y(1,4,3) - u(4,1,3) <= -0
MCLin8(4,1,3):  53 y(1,4,3) + v(4,1,3) - u(4,1,3) <= 53
MCLin9(4,1,3):  53 y(1,4,3) - u(4,1,3) >= -0
MCLin10(4,1,3):  -87 y(1,4,3) + v(4,1,3) - u(4,1,3) >= -87
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 3
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 3
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -68 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  80 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 80
MCLin4(2,3,0):  80 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -68 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -68
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -86 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  62 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 62
MCLin9(3,2,0):  62 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -86 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -86
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 1
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 2
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -74 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  74 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 74
MCLin4(2,3,1):  74 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -74 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -74
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -77 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  71 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 71
MCLin9(3,2,1):  71 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -77 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -77
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 1
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 9
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -89 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  59 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 59
MCLin4(2,3,2):  59 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -89 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -89
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -69 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  79 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 79
MCLin9(3,2,2):  79 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -69 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -69
MCLinV(2,3,3):  x(2,3) - v(2,3,3) - x(3,3) = 2
MCLinV(3,2,3):  - x(2,3) - v(3,2,3) + x(3,3) = 9
MCLin1(2,3,3):  u(2,3,3) >= -0
MCLin2(2,3,3):  -74 y(3,2,3) - u(2,3,3) <= -0
MCLin3(2,3,3):  74 y(3,2,3) + v(2,3,3) - u(2,3,3) <= 74
MCLin4(2,3,3):  74 y(3,2,3) - u(2,3,3) >= -0
MCLin5(2,3,3):  -74 y(3,2,3) + v(2,3,3) - u(2,3,3) >= -74
MCLin6(3,2,3):  u(3,2,3) >= -0
MCLin7(3,2,3):  -85 y(2,3,3) - u(3,2,3) <= -0
MCLin8(3,2,3):  63 y(2,3,3) + v(3,2,3) - u(3,2,3) <= 63
MCLin9(3,2,3):  63 y(2,3,3) - u(3,2,3) >= -0
MCLin10(3,2,3):  -85 y(2,3,3) + v(3,2,3) - u(3,2,3) >= -85
MCLinV(2,4,0):  x(2,0) - v(2,4,0) - x(4,0) = 8
MCLinV(4,2,0):  - x(2,0) - v(4,2,0) + x(4,0) = 3
MCLin1(2,4,0):  u(2,4,0) >= -0
MCLin2(2,4,0):  -59 y(4,2,0) - u(2,4,0) <= -0
MCLin3(2,4,0):  77 y(4,2,0) + v(2,4,0) - u(2,4,0) <= 77
MCLin4(2,4,0):  77 y(4,2,0) - u(2,4,0) >= -0
MCLin5(2,4,0):  -59 y(4,2,0) + v(2,4,0) - u(2,4,0) >= -59
MCLin6(4,2,0):  u(4,2,0) >= -0
MCLin7(4,2,0):  -88 y(2,4,0) - u(4,2,0) <= -0
MCLin8(4,2,0):  48 y(2,4,0) + v(4,2,0) - u(4,2,0) <= 48
MCLin9(4,2,0):  48 y(2,4,0) - u(4,2,0) >= -0
MCLin10(4,2,0):  -88 y(2,4,0) + v(4,2,0) - u(4,2,0) >= -88
MCLinV(2,4,1):  x(2,1) - v(2,4,1) - x(4,1) = 9
MCLinV(4,2,1):  - x(2,1) - v(4,2,1) + x(4,1) = 2
MCLin1(2,4,1):  u(2,4,1) >= -0
MCLin2(2,4,1):  -80 y(4,2,1) - u(2,4,1) <= -0
MCLin3(2,4,1):  56 y(4,2,1) + v(2,4,1) - u(2,4,1) <= 56
MCLin4(2,4,1):  56 y(4,2,1) - u(2,4,1) >= -0
MCLin5(2,4,1):  -80 y(4,2,1) + v(2,4,1) - u(2,4,1) >= -80
MCLin6(4,2,1):  u(4,2,1) >= -0
MCLin7(4,2,1):  -67 y(2,4,1) - u(4,2,1) <= -0
MCLin8(4,2,1):  69 y(2,4,1) + v(4,2,1) - u(4,2,1) <= 69
MCLin9(4,2,1):  69 y(2,4,1) - u(4,2,1) >= -0
MCLin10(4,2,1):  -67 y(2,4,1) + v(4,2,1) - u(4,2,1) >= -67
MCLinV(2,4,2):  x(2,2) - v(2,4,2) - x(4,2) = 1
MCLinV(4,2,2):  - x(2,2) - v(4,2,2) + x(4,2) = 9
MCLin1(2,4,2):  u(2,4,2) >= -0
MCLin2(2,4,2):  -80 y(4,2,2) - u(2,4,2) <= -0
MCLin3(2,4,2):  56 y(4,2,2) + v(2,4,2) - u(2,4,2) <= 56
MCLin4(2,4,2):  56 y(4,2,2) - u(2,4,2) >= -0
MCLin5(2,4,2):  -80 y(4,2,2) + v(2,4,2) - u(2,4,2) >= -80
MCLin6(4,2,2):  u(4,2,2) >= -0
MCLin7(4,2,2):  -66 y(2,4,2) - u(4,2,2) <= -0
MCLin8(4,2,2):  70 y(2,4,2) + v(4,2,2) - u(4,2,2) <= 70
MCLin9(4,2,2):  70 y(2,4,2) - u(4,2,2) >= -0
MCLin10(4,2,2):  -66 y(2,4,2) + v(4,2,2) - u(4,2,2) >= -66
MCLinV(2,4,3):  x(2,3) - v(2,4,3) - x(4,3) = 1
MCLinV(4,2,3):  - x(2,3) - v(4,2,3) + x(4,3) = 9
MCLin1(2,4,3):  u(2,4,3) >= -0
MCLin2(2,4,3):  -60 y(4,2,3) - u(2,4,3) <= -0
MCLin3(2,4,3):  76 y(4,2,3) + v(2,4,3) - u(2,4,3) <= 76
MCLin4(2,4,3):  76 y(4,2,3) - u(2,4,3) >= -0
MCLin5(2,4,3):  -60 y(4,2,3) + v(2,4,3) - u(2,4,3) >= -60
MCLin6(4,2,3):  u(4,2,3) >= -0
MCLin7(4,2,3):  -86 y(2,4,3) - u(4,2,3) <= -0
MCLin8(4,2,3):  50 y(2,4,3) + v(4,2,3) - u(4,2,3) <= 50
MCLin9(4,2,3):  50 y(2,4,3) - u(4,2,3) >= -0
MCLin10(4,2,3):  -86 y(2,4,3) + v(4,2,3) - u(4,2,3) >= -86
MCLinV(3,4,0):  x(3,0) - v(3,4,0) - x(4,0) = 8
MCLinV(4,3,0):  - x(3,0) - v(4,3,0) + x(4,0) = 3
MCLin1(3,4,0):  u(3,4,0) >= -0
MCLin2(3,4,0):  -76 y(4,3,0) - u(3,4,0) <= -0
MCLin3(3,4,0):  76 y(4,3,0) + v(3,4,0) - u(3,4,0) <= 76
MCLin4(3,4,0):  76 y(4,3,0) - u(3,4,0) >= -0
MCLin5(3,4,0):  -76 y(4,3,0) + v(3,4,0) - u(3,4,0) >= -76
MCLin6(4,3,0):  u(4,3,0) >= -0
MCLin7(4,3,0):  -87 y(3,4,0) - u(4,3,0) <= -0
MCLin8(4,3,0):  65 y(3,4,0) + v(4,3,0) - u(4,3,0) <= 65
MCLin9(4,3,0):  65 y(3,4,0) - u(4,3,0) >= -0
MCLin10(4,3,0):  -87 y(3,4,0) + v(4,3,0) - u(4,3,0) >= -87
MCLinV(3,4,1):  x(3,1) - v(3,4,1) - x(4,1) = 9
MCLinV(4,3,1):  - x(3,1) - v(4,3,1) + x(4,1) = 1
MCLin1(3,4,1):  u(3,4,1) >= -0
MCLin2(3,4,1):  -89 y(4,3,1) - u(3,4,1) <= -0
MCLin3(3,4,1):  63 y(4,3,1) + v(3,4,1) - u(3,4,1) <= 63
MCLin4(3,4,1):  63 y(4,3,1) - u(3,4,1) >= -0
MCLin5(3,4,1):  -89 y(4,3,1) + v(3,4,1) - u(3,4,1) >= -89
MCLin6(4,3,1):  u(4,3,1) >= -0
MCLin7(4,3,1):  -73 y(3,4,1) - u(4,3,1) <= -0
MCLin8(4,3,1):  79 y(3,4,1) + v(4,3,1) - u(4,3,1) <= 79
MCLin9(4,3,1):  79 y(3,4,1) - u(4,3,1) >= -0
MCLin10(4,3,1):  -73 y(3,4,1) + v(4,3,1) - u(4,3,1) >= -73
MCLinV(3,4,2):  x(3,2) - v(3,4,2) - x(4,2) = 1
MCLinV(4,3,2):  - x(3,2) - v(4,3,2) + x(4,2) = 1
MCLin1(3,4,2):  u(3,4,2) >= -0
MCLin2(3,4,2):  -74 y(4,3,2) - u(3,4,2) <= -0
MCLin3(3,4,2):  78 y(4,3,2) + v(3,4,2) - u(3,4,2) <= 78
MCLin4(3,4,2):  78 y(4,3,2) - u(3,4,2) >= -0
MCLin5(3,4,2):  -74 y(4,3,2) + v(3,4,2) - u(3,4,2) >= -74
MCLin6(4,3,2):  u(4,3,2) >= -0
MCLin7(4,3,2):  -80 y(3,4,2) - u(4,3,2) <= -0
MCLin8(4,3,2):  72 y(3,4,2) + v(4,3,2) - u(4,3,2) <= 72
MCLin9(4,3,2):  72 y(3,4,2) - u(4,3,2) >= -0
MCLin10(4,3,2):  -80 y(3,4,2) + v(4,3,2) - u(4,3,2) >= -80
MCLinV(3,4,3):  x(3,3) - v(3,4,3) - x(4,3) = 1
MCLinV(4,3,3):  - x(3,3) - v(4,3,3) + x(4,3) = 2
MCLin1(3,4,3):  u(3,4,3) >= -0
MCLin2(3,4,3):  -70 y(4,3,3) - u(3,4,3) <= -0
MCLin3(3,4,3):  82 y(4,3,3) + v(3,4,3) - u(3,4,3) <= 82
MCLin4(3,4,3):  82 y(4,3,3) - u(3,4,3) >= -0
MCLin5(3,4,3):  -70 y(4,3,3) + v(3,4,3) - u(3,4,3) >= -70
MCLin6(4,3,3):  u(4,3,3) >= -0
MCLin7(4,3,3):  -85 y(3,4,3) - u(4,3,3) <= -0
MCLin8(4,3,3):  67 y(3,4,3) + v(4,3,3) - u(4,3,3) <= 67
MCLin9(4,3,3):  67 y(3,4,3) - u(4,3,3) >= -0
MCLin10(4,3,3):  -85 y(3,4,3) + v(4,3,3) - u(4,3,3) >= -85
makespan(0):  C - x(0,3) >= 1
makespan(1):  C - x(1,2) >= 2
makespan(2):  C - x(2,0) >= 3
makespan(3):  C - x(3,2) >= 1
makespan(4):  C - x(4,1) >= 9
cut_clique_0(01234,0):  0.04044 x(1,0) + 0.04412 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(0134,0):  0.04044 x(1,0) + 0.04412 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(1234,0):  0.04044 x(1,0) + 0.04412 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(0234,0):  0.09804 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(0123,0):  0.04167 x(1,0) + 0.25000 x(3,0) >= 1
cut_clique_0(0124,0):  0.08333 x(1,0) + 0.25000 x(4,0) >= 1
cut_clique_0(134,0):  0.04044 x(1,0) + 0.04412 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(034,0):  0.09804 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(013,0):  0.04167 x(1,0) + 0.25000 x(3,0) >= 1
cut_clique_0(014,0):  0.08333 x(1,0) + 0.25000 x(4,0) >= 1
cut_clique_0(234,0):  0.09804 x(3,0) + 0.11765 x(4,0) >= 1
cut_clique_0(123,0):  0.04167 x(1,0) + 0.25000 x(3,0) >= 1
cut_clique_0(124,0):  0.08333 x(1,0) + 0.25000 x(4,0) >= 1
cut_clique_0(01234,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(0123,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(0134,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(0124,1):  0.02857 x(0,1) + 0.20000 x(1,1) >= 1
cut_clique_0(1234,1):  0.00840 x(2,1) + 0.08403 x(4,1) >= 1
cut_clique_0(0234,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(013,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(124,1):  0.00840 x(2,1) + 0.08403 x(4,1) >= 1
cut_clique_0(012,1):  0.02857 x(0,1) + 0.20000 x(1,1) >= 1
cut_clique_0(023,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(014,1):  0.02857 x(0,1) + 0.20000 x(1,1) >= 1
cut_clique_0(024,1):  0.00840 x(2,1) + 0.08403 x(4,1) >= 1
cut_clique_0(034,1):  x(0,1) + 0.20000 x(3,1) >= 1
cut_clique_0(234,1):  0.00840 x(2,1) + 0.08403 x(4,1) >= 1
cut_clique_0(01234,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(3,2) >= 1
cut_clique_0(0234,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(3,2) >= 1
cut_clique_0(0134,2):  0.05512 x(0,2) + 0.02520 x(3,2) + 0.03307 x(4,2) >= 1
cut_clique_0(0123,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(3,2) >= 1
cut_clique_0(1234,2):  0.04762 x(2,2) + 0.11111 x(3,2) >= 1
cut_clique_0(0124,2):  0.11111 x(0,2) + 0.03704 x(2,2) >= 1
cut_clique_0(034,2):  0.05512 x(0,2) + 0.02520 x(3,2) + 0.03307 x(4,2) >= 1
cut_clique_0(023,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(3,2) >= 1
cut_clique_0(024,2):  0.11111 x(0,2) + 0.03704 x(2,2) >= 1
cut_clique_0(234,2):  0.04762 x(2,2) + 0.11111 x(3,2) >= 1
cut_clique_0(013,2):  0.11111 x(0,2) + 0.03704 x(3,2) >= 1
cut_clique_0(014,2):  0.04000 x(0,2) + 0.06667 x(4,2) >= 1
cut_clique_0(012,2):  0.11111 x(0,2) + 0.03704 x(2,2) >= 1
cut_clique_0(123,2):  0.04762 x(2,2) + 0.11111 x(3,2) >= 1
cut_clique_0(01234,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(1234,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(0123,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(0124,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(123,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(124,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_0(012,3):  0.02281 x(1,3) + 0.04943 x(2,3) >= 1
cut_clique_1(01234,0):  0.01236 x(0,0) + 0.01236 x(1,0) + 0.02163 x(2,0) + 0.00927 x(3,0) + 0.05149 x(4,0) >= 1
cut_clique_1(0134,0):  0.03659 x(0,0) + 0.01829 x(1,0) + 0.01372 x(3,0) + 0.07165 x(4,0) >= 1
cut_clique_1(1234,0):  0.04545 x(1,0) + 0.03409 x(3,0) + 0.14773 x(4,0) >= 1
cut_clique_1(134,0):  0.04545 x(1,0) + 0.03409 x(3,0) + 0.14773 x(4,0) >= 1
cut_clique_1(01234,2):  0.07292 x(0,2) + 0.02679 x(2,2) + 0.01042 x(3,2) + 0.01042 x(4,2) >= 1
cut_clique_1(0234,2):  0.07292 x(0,2) + 0.02679 x(2,2) + 0.01042 x(3,2) + 0.01042 x(4,2) >= 1
cut_clique_1(0124,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(4,2) >= 1
cut_clique_1(024,2):  0.08861 x(0,2) + 0.03116 x(2,2) + 0.01266 x(4,2) >= 1
cut_clique_1(01234,3):  0.00210 x(0,3) + 0.02102 x(1,3) + 0.04695 x(2,3) >= 1
cut_clique_1(0123,3):  0.00210 x(0,3) + 0.02102 x(1,3) + 0.04695 x(2,3) >= 1
cut_clique_1(0124,3):  0.00210 x(0,3) + 0.02102 x(1,3) + 0.04695 x(2,3) >= 1
cut_clique_1(012,3):  0.00210 x(0,3) + 0.02102 x(1,3) + 0.04695 x(2,3) >= 1
cut_clique_0:  0.00719 x(2,1) + 0.07194 x(4,1) + 0.01439 x(0,2) >= 1
cut_clique_1:  0.07813 x(4,0) + 0.00078 x(0,1) + 0.01563 x(1,1) + 0.01042 x(2,1) + 0.04687 x(4,1) >= 1
Bounds
 0 <= C <= 89
 12 <= x(0,0) <= 80
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -62 <= v(0,1,0) <= 76
 -88 <= v(1,0,0) <= 50
 -890 <= u(0,1,0)
 -890 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -77 <= v(0,2,0) <= 57
 -68 <= v(2,0,0) <= 66
 -890 <= u(0,2,0)
 -890 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -76 <= v(0,3,0) <= 74
 -85 <= v(3,0,0) <= 65
 -890 <= u(0,3,0)
 -890 <= u(3,0,0)
 0 <= y(0,4,0) <= 1
 0 <= y(4,0,0) <= 1
 -67 <= v(0,4,0) <= 71
 -87 <= v(4,0,0) <= 51
 -890 <= u(0,4,0)
 -890 <= u(4,0,0)
 0 <= x(1,0) <= 70
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -89 <= v(1,2,0) <= 47
 -54 <= v(2,1,0) <= 82
 -890 <= u(1,2,0)
 -890 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -88 <= v(1,3,0) <= 64
 -71 <= v(3,1,0) <= 81
 -890 <= u(1,3,0)
 -890 <= u(3,1,0)
 0 <= y(1,4,0) <= 1
 0 <= y(4,1,0) <= 1
 -79 <= v(1,4,0) <= 61
 -73 <= v(4,1,0) <= 67
 -890 <= u(1,4,0)
 -890 <= u(4,1,0)
 20 <= x(2,0) <= 86
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -68 <= v(2,3,0) <= 80
 -86 <= v(3,2,0) <= 62
 -890 <= u(2,3,0)
 -890 <= u(3,2,0)
 0 <= y(2,4,0) <= 1
 0 <= y(4,2,0) <= 1
 -59 <= v(2,4,0) <= 77
 -88 <= v(4,2,0) <= 48
 -890 <= u(2,4,0)
 -890 <= u(4,2,0)
 3 <= x(3,0) <= 85
 0 <= y(3,4,0) <= 1
 0 <= y(4,3,0) <= 1
 -76 <= v(3,4,0) <= 76
 -87 <= v(4,3,0) <= 65
 -890 <= u(3,4,0)
 -890 <= u(4,3,0)
 1 <= x(4,0) <= 71
 0 <= x(0,1) <= 68
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -77 <= v(0,1,1) <= 61
 -69 <= v(1,0,1) <= 69
 -890 <= u(0,1,1)
 -890 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -77 <= v(0,2,1) <= 57
 -64 <= v(2,0,1) <= 70
 -890 <= u(0,2,1)
 -890 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -83 <= v(0,3,1) <= 67
 -73 <= v(3,0,1) <= 77
 -890 <= u(0,3,1)
 -890 <= u(3,0,1)
 0 <= y(0,4,1) <= 1
 0 <= y(4,0,1) <= 1
 -89 <= v(0,4,1) <= 49
 -63 <= v(4,0,1) <= 75
 -890 <= u(0,4,1)
 -890 <= u(4,0,1)
 4 <= x(1,1) <= 74
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -73 <= v(1,2,1) <= 63
 -68 <= v(2,1,1) <= 68
 -890 <= u(1,2,1)
 -890 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -79 <= v(1,3,1) <= 73
 -77 <= v(3,1,1) <= 75
 -890 <= u(1,3,1)
 -890 <= u(3,1,1)
 0 <= y(1,4,1) <= 1
 0 <= y(4,1,1) <= 1
 -85 <= v(1,4,1) <= 55
 -67 <= v(4,1,1) <= 73
 -890 <= u(1,4,1)
 -890 <= u(4,1,1)
 9 <= x(2,1) <= 75
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -74 <= v(2,3,1) <= 74
 -77 <= v(3,2,1) <= 71
 -890 <= u(2,3,1)
 -890 <= u(3,2,1)
 0 <= y(2,4,1) <= 1
 0 <= y(4,2,1) <= 1
 -80 <= v(2,4,1) <= 56
 -67 <= v(4,2,1) <= 69
 -890 <= u(2,4,1)
 -890 <= u(4,2,1)
 0 <= x(3,1) <= 82
 0 <= y(3,4,1) <= 1
 0 <= y(4,3,1) <= 1
 -89 <= v(3,4,1) <= 63
 -73 <= v(4,3,1) <= 79
 -890 <= u(3,4,1)
 -890 <= u(4,3,1)
 10 <= x(4,1) <= 80
 5 <= x(0,2) <= 73
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -84 <= v(0,1,2) <= 54
 -63 <= v(1,0,2) <= 75
 -890 <= u(0,1,2)
 -890 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -70 <= v(0,2,2) <= 64
 -80 <= v(2,0,2) <= 54
 -890 <= u(0,2,2)
 -890 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -84 <= v(0,3,2) <= 66
 -74 <= v(3,0,2) <= 76
 -890 <= u(0,3,2)
 -890 <= u(3,0,2)
 0 <= y(0,4,2) <= 1
 0 <= y(4,0,2) <= 1
 -75 <= v(0,4,2) <= 63
 -71 <= v(4,0,2) <= 67
 -890 <= u(0,4,2)
 -890 <= u(4,0,2)
 17 <= x(1,2) <= 87
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -58 <= v(1,2,2) <= 78
 -89 <= v(2,1,2) <= 47
 -890 <= u(1,2,2)
 -890 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -72 <= v(1,3,2) <= 80
 -83 <= v(3,1,2) <= 69
 -890 <= u(1,3,2)
 -890 <= u(3,1,2)
 0 <= y(1,4,2) <= 1
 0 <= y(4,1,2) <= 1
 -63 <= v(1,4,2) <= 77
 -80 <= v(4,1,2) <= 60
 -890 <= u(1,4,2)
 -890 <= u(4,1,2)
 0 <= x(2,2) <= 66
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -89 <= v(2,3,2) <= 59
 -69 <= v(3,2,2) <= 79
 -890 <= u(2,3,2)
 -890 <= u(3,2,2)
 0 <= y(2,4,2) <= 1
 0 <= y(4,2,2) <= 1
 -80 <= v(2,4,2) <= 56
 -66 <= v(4,2,2) <= 70
 -890 <= u(2,4,2)
 -890 <= u(4,2,2)
 6 <= x(3,2) <= 88
 0 <= y(3,4,2) <= 1
 0 <= y(4,3,2) <= 1
 -74 <= v(3,4,2) <= 78
 -80 <= v(4,3,2) <= 72
 -890 <= u(3,4,2)
 -890 <= u(4,3,2)
 9 <= x(4,2) <= 79
 20 <= x(0,3) <= 88
 0 <= y(0,1,3) <= 1
 0 <= y(1,0,3) <= 1
 -67 <= v(0,1,3) <= 71
 -82 <= v(1,0,3) <= 56
 -890 <= u(0,1,3)
 -890 <= u(1,0,3)
 0 <= y(0,2,3) <= 1
 0 <= y(2,0,3) <= 1
 -66 <= v(0,2,3) <= 68
 -78 <= v(2,0,3) <= 56
 -890 <= u(0,2,3)
 -890 <= u(2,0,3)
 0 <= y(0,3,3) <= 1
 0 <= y(3,0,3) <= 1
 -65 <= v(0,3,3) <= 85
 -88 <= v(3,0,3) <= 62
 -890 <= u(0,3,3)
 -890 <= u(3,0,3)
 0 <= y(0,4,3) <= 1
 0 <= y(4,0,3) <= 1
 -51 <= v(0,4,3) <= 87
 -89 <= v(4,0,3) <= 49
 -890 <= u(0,4,3)
 -890 <= u(4,0,3)
 7 <= x(1,3) <= 77
 0 <= y(1,2,3) <= 1
 0 <= y(2,1,3) <= 1
 -79 <= v(1,2,3) <= 57
 -76 <= v(2,1,3) <= 60
 -890 <= u(1,2,3)
 -890 <= u(2,1,3)
 0 <= y(1,3,3) <= 1
 0 <= y(3,1,3) <= 1
 -78 <= v(1,3,3) <= 74
 -86 <= v(3,1,3) <= 66
 -890 <= u(1,3,3)
 -890 <= u(3,1,3)
 0 <= y(1,4,3) <= 1
 0 <= y(4,1,3) <= 1
 -64 <= v(1,4,3) <= 76
 -87 <= v(4,1,3) <= 53
 -890 <= u(1,4,3)
 -890 <= u(4,1,3)
 11 <= x(2,3) <= 77
 0 <= y(2,3,3) <= 1
 0 <= y(3,2,3) <= 1
 -74 <= v(2,3,3) <= 74
 -85 <= v(3,2,3) <= 63
 -890 <= u(2,3,3)
 -890 <= u(3,2,3)
 0 <= y(2,4,3) <= 1
 0 <= y(4,2,3) <= 1
 -60 <= v(2,4,3) <= 76
 -86 <= v(4,2,3) <= 50
 -890 <= u(2,4,3)
 -890 <= u(4,2,3)
 1 <= x(3,3) <= 83
 0 <= y(3,4,3) <= 1
 0 <= y(4,3,3) <= 1
 -70 <= v(3,4,3) <= 82
 -85 <= v(4,3,3) <= 67
 -890 <= u(3,4,3)
 -890 <= u(4,3,3)
 0 <= x(4,3) <= 70
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
