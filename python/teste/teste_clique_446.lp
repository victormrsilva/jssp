\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,3) >= 9
ord(0,2):  x(0,2) - x(0,3) >= 8
ord(0,3):  x(0,0) - x(0,2) >= 9
ord(1,1):  - x(1,2) + x(1,3) >= 7
ord(1,2):  x(1,0) - x(1,3) >= 9
ord(1,3):  - x(1,0) + x(1,1) >= 9
ord(2,1):  x(2,1) - x(2,3) >= 10
ord(2,2):  x(2,0) - x(2,1) >= 2
ord(2,3):  - x(2,0) + x(2,2) >= 10
ord(3,1):  x(3,1) - x(3,3) >= 2
ord(3,2):  x(3,0) - x(3,1) >= 10
ord(3,3):  - x(3,0) + x(3,2) >= 7
ord(4,1):  - x(4,1) + x(4,2) >= 10
ord(4,2):  - x(4,2) + x(4,3) >= 2
ord(4,3):  x(4,0) - x(4,3) >= 6
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
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 9
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 6
MCLin1(0,1,0):  u(0,1,0) >= -0
MCLin2(0,1,0):  -96 y(1,0,0) - u(0,1,0) <= -0
MCLin3(0,1,0):  93 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 93
MCLin4(0,1,0):  93 y(1,0,0) - u(0,1,0) >= -0
MCLin5(0,1,0):  -96 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -96
MCLin6(1,0,0):  u(1,0,0) >= -0
MCLin7(1,0,0):  -108 y(0,1,0) - u(1,0,0) <= -0
MCLin8(1,0,0):  81 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 81
MCLin9(1,0,0):  81 y(0,1,0) - u(1,0,0) >= -0
MCLin10(1,0,0):  -108 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -108
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 2
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 9
MCLin1(0,1,1):  u(0,1,1) >= -0
MCLin2(0,1,1):  -124 y(1,0,1) - u(0,1,1) <= -0
MCLin3(0,1,1):  65 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 65
MCLin4(0,1,1):  65 y(1,0,1) - u(0,1,1) >= -0
MCLin5(0,1,1):  -124 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -124
MCLin6(1,0,1):  u(1,0,1) >= -0
MCLin7(1,0,1):  -76 y(0,1,1) - u(1,0,1) <= -0
MCLin8(1,0,1):  113 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 113
MCLin9(1,0,1):  113 y(0,1,1) - u(1,0,1) >= -0
MCLin10(1,0,1):  -76 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -76
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 7
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 9
MCLin1(0,1,2):  u(0,1,2) >= -0
MCLin2(0,1,2):  -87 y(1,0,2) - u(0,1,2) <= -0
MCLin3(0,1,2):  102 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 102
MCLin4(0,1,2):  102 y(1,0,2) - u(0,1,2) >= -0
MCLin5(0,1,2):  -87 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -87
MCLin6(1,0,2):  u(1,0,2) >= -0
MCLin7(1,0,2):  -118 y(0,1,2) - u(1,0,2) <= -0
MCLin8(1,0,2):  71 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 71
MCLin9(1,0,2):  71 y(0,1,2) - u(1,0,2) >= -0
MCLin10(1,0,2):  -118 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -118
MCLinV(0,1,3):  x(0,3) - v(0,1,3) - x(1,3) = 9
MCLinV(1,0,3):  - x(0,3) - v(1,0,3) + x(1,3) = 8
MCLin1(0,1,3):  u(0,1,3) >= -0
MCLin2(0,1,3):  -104 y(1,0,3) - u(0,1,3) <= -0
MCLin3(0,1,3):  85 y(1,0,3) + v(0,1,3) - u(0,1,3) <= 85
MCLin4(0,1,3):  85 y(1,0,3) - u(0,1,3) >= -0
MCLin5(0,1,3):  -104 y(1,0,3) + v(0,1,3) - u(0,1,3) >= -104
MCLin6(1,0,3):  u(1,0,3) >= -0
MCLin7(1,0,3):  -102 y(0,1,3) - u(1,0,3) <= -0
MCLin8(1,0,3):  87 y(0,1,3) + v(1,0,3) - u(1,0,3) <= 87
MCLin9(1,0,3):  87 y(0,1,3) - u(1,0,3) >= -0
MCLin10(1,0,3):  -102 y(0,1,3) + v(1,0,3) - u(1,0,3) >= -102
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 10
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 6
MCLin1(0,2,0):  u(0,2,0) >= -0
MCLin2(0,2,0):  -97 y(2,0,0) - u(0,2,0) <= -0
MCLin3(0,2,0):  96 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 96
MCLin4(0,2,0):  96 y(2,0,0) - u(0,2,0) >= -0
MCLin5(0,2,0):  -97 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -97
MCLin6(2,0,0):  u(2,0,0) >= -0
MCLin7(2,0,0):  -112 y(0,2,0) - u(2,0,0) <= -0
MCLin8(2,0,0):  81 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 81
MCLin9(2,0,0):  81 y(0,2,0) - u(2,0,0) >= -0
MCLin10(2,0,0):  -112 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -112
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 2
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 9
MCLin1(0,2,1):  u(0,2,1) >= -0
MCLin2(0,2,1):  -113 y(2,0,1) - u(0,2,1) <= -0
MCLin3(0,2,1):  80 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 80
MCLin4(0,2,1):  80 y(2,0,1) - u(0,2,1) >= -0
MCLin5(0,2,1):  -113 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -113
MCLin6(2,0,1):  u(2,0,1) >= -0
MCLin7(2,0,1):  -91 y(0,2,1) - u(2,0,1) <= -0
MCLin8(2,0,1):  102 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 102
MCLin9(2,0,1):  102 y(0,2,1) - u(2,0,1) >= -0
MCLin10(2,0,1):  -91 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -91
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 1
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 9
MCLin1(0,2,2):  u(0,2,2) >= -0
MCLin2(0,2,2):  -107 y(2,0,2) - u(0,2,2) <= -0
MCLin3(0,2,2):  86 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 86
MCLin4(0,2,2):  86 y(2,0,2) - u(0,2,2) >= -0
MCLin5(0,2,2):  -107 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -107
MCLin6(2,0,2):  u(2,0,2) >= -0
MCLin7(2,0,2):  -96 y(0,2,2) - u(2,0,2) <= -0
MCLin8(2,0,2):  97 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 97
MCLin9(2,0,2):  97 y(0,2,2) - u(2,0,2) >= -0
MCLin10(2,0,2):  -96 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -96
MCLinV(0,2,3):  x(0,3) - v(0,2,3) - x(2,3) = 10
MCLinV(2,0,3):  - x(0,3) - v(2,0,3) + x(2,3) = 8
MCLin1(0,2,3):  u(0,2,3) >= -0
MCLin2(0,2,3):  -102 y(2,0,3) - u(0,2,3) <= -0
MCLin3(0,2,3):  91 y(2,0,3) + v(0,2,3) - u(0,2,3) <= 91
MCLin4(0,2,3):  91 y(2,0,3) - u(0,2,3) >= -0
MCLin5(0,2,3):  -102 y(2,0,3) + v(0,2,3) - u(0,2,3) >= -102
MCLin6(2,0,3):  u(2,0,3) >= -0
MCLin7(2,0,3):  -109 y(0,2,3) - u(2,0,3) <= -0
MCLin8(2,0,3):  84 y(0,2,3) + v(2,0,3) - u(2,0,3) <= 84
MCLin9(2,0,3):  84 y(0,2,3) - u(2,0,3) >= -0
MCLin10(2,0,3):  -109 y(0,2,3) + v(2,0,3) - u(2,0,3) >= -109
MCLinV(0,3,0):  x(0,0) - v(0,3,0) - x(3,0) = 7
MCLinV(3,0,0):  - x(0,0) - v(3,0,0) + x(3,0) = 6
MCLin1(0,3,0):  u(0,3,0) >= -0
MCLin2(0,3,0):  -97 y(3,0,0) - u(0,3,0) <= -0
MCLin3(0,3,0):  99 y(3,0,0) + v(0,3,0) - u(0,3,0) <= 99
MCLin4(0,3,0):  99 y(3,0,0) - u(0,3,0) >= -0
MCLin5(0,3,0):  -97 y(3,0,0) + v(0,3,0) - u(0,3,0) >= -97
MCLin6(3,0,0):  u(3,0,0) >= -0
MCLin7(3,0,0):  -112 y(0,3,0) - u(3,0,0) <= -0
MCLin8(3,0,0):  84 y(0,3,0) + v(3,0,0) - u(3,0,0) <= 84
MCLin9(3,0,0):  84 y(0,3,0) - u(3,0,0) >= -0
MCLin10(3,0,0):  -112 y(0,3,0) + v(3,0,0) - u(3,0,0) >= -112
MCLinV(0,3,1):  x(0,1) - v(0,3,1) - x(3,1) = 10
MCLinV(3,0,1):  - x(0,1) - v(3,0,1) + x(3,1) = 9
MCLin1(0,3,1):  u(0,3,1) >= -0
MCLin2(0,3,1):  -116 y(3,0,1) - u(0,3,1) <= -0
MCLin3(0,3,1):  80 y(3,0,1) + v(0,3,1) - u(0,3,1) <= 80
MCLin4(0,3,1):  80 y(3,0,1) - u(0,3,1) >= -0
MCLin5(0,3,1):  -116 y(3,0,1) + v(0,3,1) - u(0,3,1) >= -116
MCLin6(3,0,1):  u(3,0,1) >= -0
MCLin7(3,0,1):  -99 y(0,3,1) - u(3,0,1) <= -0
MCLin8(3,0,1):  97 y(0,3,1) + v(3,0,1) - u(3,0,1) <= 97
MCLin9(3,0,1):  97 y(0,3,1) - u(3,0,1) >= -0
MCLin10(3,0,1):  -99 y(0,3,1) + v(3,0,1) - u(3,0,1) >= -99
MCLinV(0,3,2):  x(0,2) - v(0,3,2) - x(3,2) = 1
MCLinV(3,0,2):  - x(0,2) - v(3,0,2) + x(3,2) = 9
MCLin1(0,3,2):  u(0,3,2) >= -0
MCLin2(0,3,2):  -107 y(3,0,2) - u(0,3,2) <= -0
MCLin3(0,3,2):  89 y(3,0,2) + v(0,3,2) - u(0,3,2) <= 89
MCLin4(0,3,2):  89 y(3,0,2) - u(0,3,2) >= -0
MCLin5(0,3,2):  -107 y(3,0,2) + v(0,3,2) - u(0,3,2) >= -107
MCLin6(3,0,2):  u(3,0,2) >= -0
MCLin7(3,0,2):  -99 y(0,3,2) - u(3,0,2) <= -0
MCLin8(3,0,2):  97 y(0,3,2) + v(3,0,2) - u(3,0,2) <= 97
MCLin9(3,0,2):  97 y(0,3,2) - u(3,0,2) >= -0
MCLin10(3,0,2):  -99 y(0,3,2) + v(3,0,2) - u(3,0,2) >= -99
MCLinV(0,3,3):  x(0,3) - v(0,3,3) - x(3,3) = 2
MCLinV(3,0,3):  - x(0,3) - v(3,0,3) + x(3,3) = 8
MCLin1(0,3,3):  u(0,3,3) >= -0
MCLin2(0,3,3):  -97 y(3,0,3) - u(0,3,3) <= -0
MCLin3(0,3,3):  99 y(3,0,3) + v(0,3,3) - u(0,3,3) <= 99
MCLin4(0,3,3):  99 y(3,0,3) - u(0,3,3) >= -0
MCLin5(0,3,3):  -97 y(3,0,3) + v(0,3,3) - u(0,3,3) >= -97
MCLin6(3,0,3):  u(3,0,3) >= -0
MCLin7(3,0,3):  -109 y(0,3,3) - u(3,0,3) <= -0
MCLin8(3,0,3):  87 y(0,3,3) + v(3,0,3) - u(3,0,3) <= 87
MCLin9(3,0,3):  87 y(0,3,3) - u(3,0,3) >= -0
MCLin10(3,0,3):  -109 y(0,3,3) + v(3,0,3) - u(3,0,3) >= -109
MCLinV(0,4,0):  x(0,0) - v(0,4,0) - x(4,0) = 4
MCLinV(4,0,0):  - x(0,0) - v(4,0,0) + x(4,0) = 6
MCLin1(0,4,0):  u(0,4,0) >= -0
MCLin2(0,4,0):  -98 y(4,0,0) - u(0,4,0) <= -0
MCLin3(0,4,0):  96 y(4,0,0) + v(0,4,0) - u(0,4,0) <= 96
MCLin4(0,4,0):  96 y(4,0,0) - u(0,4,0) >= -0
MCLin5(0,4,0):  -98 y(4,0,0) + v(0,4,0) - u(0,4,0) >= -98
MCLin6(4,0,0):  u(4,0,0) >= -0
MCLin7(4,0,0):  -106 y(0,4,0) - u(4,0,0) <= -0
MCLin8(4,0,0):  88 y(0,4,0) + v(4,0,0) - u(4,0,0) <= 88
MCLin9(4,0,0):  88 y(0,4,0) - u(4,0,0) >= -0
MCLin10(4,0,0):  -106 y(0,4,0) + v(4,0,0) - u(4,0,0) >= -106
MCLinV(0,4,1):  x(0,1) - v(0,4,1) - x(4,1) = 10
MCLinV(4,0,1):  - x(0,1) - v(4,0,1) + x(4,1) = 9
MCLin1(0,4,1):  u(0,4,1) >= -0
MCLin2(0,4,1):  -112 y(4,0,1) - u(0,4,1) <= -0
MCLin3(0,4,1):  82 y(4,0,1) + v(0,4,1) - u(0,4,1) <= 82
MCLin4(0,4,1):  82 y(4,0,1) - u(0,4,1) >= -0
MCLin5(0,4,1):  -112 y(4,0,1) + v(0,4,1) - u(0,4,1) >= -112
MCLin6(4,0,1):  u(4,0,1) >= -0
MCLin7(4,0,1):  -101 y(0,4,1) - u(4,0,1) <= -0
MCLin8(4,0,1):  93 y(0,4,1) + v(4,0,1) - u(4,0,1) <= 93
MCLin9(4,0,1):  93 y(0,4,1) - u(4,0,1) >= -0
MCLin10(4,0,1):  -101 y(0,4,1) + v(4,0,1) - u(4,0,1) >= -101
MCLinV(0,4,2):  x(0,2) - v(0,4,2) - x(4,2) = 2
MCLinV(4,0,2):  - x(0,2) - v(4,0,2) + x(4,2) = 9
MCLin1(0,4,2):  u(0,4,2) >= -0
MCLin2(0,4,2):  -97 y(4,0,2) - u(0,4,2) <= -0
MCLin3(0,4,2):  97 y(4,0,2) + v(0,4,2) - u(0,4,2) <= 97
MCLin4(0,4,2):  97 y(4,0,2) - u(0,4,2) >= -0
MCLin5(0,4,2):  -97 y(4,0,2) + v(0,4,2) - u(0,4,2) >= -97
MCLin6(4,0,2):  u(4,0,2) >= -0
MCLin7(4,0,2):  -108 y(0,4,2) - u(4,0,2) <= -0
MCLin8(4,0,2):  86 y(0,4,2) + v(4,0,2) - u(4,0,2) <= 86
MCLin9(4,0,2):  86 y(0,4,2) - u(4,0,2) >= -0
MCLin10(4,0,2):  -108 y(0,4,2) + v(4,0,2) - u(4,0,2) >= -108
MCLinV(0,4,3):  x(0,3) - v(0,4,3) - x(4,3) = 6
MCLinV(4,0,3):  - x(0,3) - v(4,0,3) + x(4,3) = 8
MCLin1(0,4,3):  u(0,4,3) >= -0
MCLin2(0,4,3):  -111 y(4,0,3) - u(0,4,3) <= -0
MCLin3(0,4,3):  83 y(4,0,3) + v(0,4,3) - u(0,4,3) <= 83
MCLin4(0,4,3):  83 y(4,0,3) - u(0,4,3) >= -0
MCLin5(0,4,3):  -111 y(4,0,3) + v(0,4,3) - u(0,4,3) >= -111
MCLin6(4,0,3):  u(4,0,3) >= -0
MCLin7(4,0,3):  -97 y(0,4,3) - u(4,0,3) <= -0
MCLin8(4,0,3):  97 y(0,4,3) + v(4,0,3) - u(4,0,3) <= 97
MCLin9(4,0,3):  97 y(0,4,3) - u(4,0,3) >= -0
MCLin10(4,0,3):  -97 y(0,4,3) + v(4,0,3) - u(4,0,3) >= -97
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 10
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 9
MCLin1(1,2,0):  u(1,2,0) >= -0
MCLin2(1,2,0):  -107 y(2,1,0) - u(1,2,0) <= -0
MCLin3(1,2,0):  91 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 91
MCLin4(1,2,0):  91 y(2,1,0) - u(1,2,0) >= -0
MCLin5(1,2,0):  -107 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -107
MCLin6(2,1,0):  u(2,1,0) >= -0
MCLin7(2,1,0):  -110 y(1,2,0) - u(2,1,0) <= -0
MCLin8(2,1,0):  88 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 88
MCLin9(2,1,0):  88 y(1,2,0) - u(2,1,0) >= -0
MCLin10(2,1,0):  -110 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -110
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 2
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 2
MCLin1(1,2,1):  u(1,2,1) >= -0
MCLin2(1,2,1):  -88 y(2,1,1) - u(1,2,1) <= -0
MCLin3(1,2,1):  110 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 110
MCLin4(1,2,1):  110 y(2,1,1) - u(1,2,1) >= -0
MCLin5(1,2,1):  -88 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -88
MCLin6(2,1,1):  u(2,1,1) >= -0
MCLin7(2,1,1):  -114 y(1,2,1) - u(2,1,1) <= -0
MCLin8(2,1,1):  84 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 84
MCLin9(2,1,1):  84 y(1,2,1) - u(2,1,1) >= -0
MCLin10(2,1,1):  -114 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -114
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 1
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 7
MCLin1(1,2,2):  u(1,2,2) >= -0
MCLin2(1,2,2):  -124 y(2,1,2) - u(1,2,2) <= -0
MCLin3(1,2,2):  74 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 74
MCLin4(1,2,2):  74 y(2,1,2) - u(1,2,2) >= -0
MCLin5(1,2,2):  -124 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -124
MCLin6(2,1,2):  u(2,1,2) >= -0
MCLin7(2,1,2):  -82 y(1,2,2) - u(2,1,2) <= -0
MCLin8(2,1,2):  116 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 116
MCLin9(2,1,2):  116 y(1,2,2) - u(2,1,2) >= -0
MCLin10(2,1,2):  -82 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -82
MCLinV(1,2,3):  x(1,3) - v(1,2,3) - x(2,3) = 10
MCLinV(2,1,3):  - x(1,3) - v(2,1,3) + x(2,3) = 9
MCLin1(1,2,3):  u(1,2,3) >= -0
MCLin2(1,2,3):  -104 y(2,1,3) - u(1,2,3) <= -0
MCLin3(1,2,3):  94 y(2,1,3) + v(1,2,3) - u(1,2,3) <= 94
MCLin4(1,2,3):  94 y(2,1,3) - u(1,2,3) >= -0
MCLin5(1,2,3):  -104 y(2,1,3) + v(1,2,3) - u(1,2,3) >= -104
MCLin6(2,1,3):  u(2,1,3) >= -0
MCLin7(2,1,3):  -113 y(1,2,3) - u(2,1,3) <= -0
MCLin8(2,1,3):  85 y(1,2,3) + v(2,1,3) - u(2,1,3) <= 85
MCLin9(2,1,3):  85 y(1,2,3) - u(2,1,3) >= -0
MCLin10(2,1,3):  -113 y(1,2,3) + v(2,1,3) - u(2,1,3) >= -113
MCLinV(1,3,0):  x(1,0) - v(1,3,0) - x(3,0) = 7
MCLinV(3,1,0):  - x(1,0) - v(3,1,0) + x(3,0) = 9
MCLin1(1,3,0):  u(1,3,0) >= -0
MCLin2(1,3,0):  -107 y(3,1,0) - u(1,3,0) <= -0
MCLin3(1,3,0):  94 y(3,1,0) + v(1,3,0) - u(1,3,0) <= 94
MCLin4(1,3,0):  94 y(3,1,0) - u(1,3,0) >= -0
MCLin5(1,3,0):  -107 y(3,1,0) + v(1,3,0) - u(1,3,0) >= -107
MCLin6(3,1,0):  u(3,1,0) >= -0
MCLin7(3,1,0):  -110 y(1,3,0) - u(3,1,0) <= -0
MCLin8(3,1,0):  91 y(1,3,0) + v(3,1,0) - u(3,1,0) <= 91
MCLin9(3,1,0):  91 y(1,3,0) - u(3,1,0) >= -0
MCLin10(3,1,0):  -110 y(1,3,0) + v(3,1,0) - u(3,1,0) >= -110
MCLinV(1,3,1):  x(1,1) - v(1,3,1) - x(3,1) = 10
MCLinV(3,1,1):  - x(1,1) - v(3,1,1) + x(3,1) = 2
MCLin1(1,3,1):  u(1,3,1) >= -0
MCLin2(1,3,1):  -91 y(3,1,1) - u(1,3,1) <= -0
MCLin3(1,3,1):  110 y(3,1,1) + v(1,3,1) - u(1,3,1) <= 110
MCLin4(1,3,1):  110 y(3,1,1) - u(1,3,1) >= -0
MCLin5(1,3,1):  -91 y(3,1,1) + v(1,3,1) - u(1,3,1) >= -91
MCLin6(3,1,1):  u(3,1,1) >= -0
MCLin7(3,1,1):  -122 y(1,3,1) - u(3,1,1) <= -0
MCLin8(3,1,1):  79 y(1,3,1) + v(3,1,1) - u(3,1,1) <= 79
MCLin9(3,1,1):  79 y(1,3,1) - u(3,1,1) >= -0
MCLin10(3,1,1):  -122 y(1,3,1) + v(3,1,1) - u(3,1,1) >= -122
MCLinV(1,3,2):  x(1,2) - v(1,3,2) - x(3,2) = 1
MCLinV(3,1,2):  - x(1,2) - v(3,1,2) + x(3,2) = 7
MCLin1(1,3,2):  u(1,3,2) >= -0
MCLin2(1,3,2):  -124 y(3,1,2) - u(1,3,2) <= -0
MCLin3(1,3,2):  77 y(3,1,2) + v(1,3,2) - u(1,3,2) <= 77
MCLin4(1,3,2):  77 y(3,1,2) - u(1,3,2) >= -0
MCLin5(1,3,2):  -124 y(3,1,2) + v(1,3,2) - u(1,3,2) >= -124
MCLin6(3,1,2):  u(3,1,2) >= -0
MCLin7(3,1,2):  -85 y(1,3,2) - u(3,1,2) <= -0
MCLin8(3,1,2):  116 y(1,3,2) + v(3,1,2) - u(3,1,2) <= 116
MCLin9(3,1,2):  116 y(1,3,2) - u(3,1,2) >= -0
MCLin10(3,1,2):  -85 y(1,3,2) + v(3,1,2) - u(3,1,2) >= -85
MCLinV(1,3,3):  x(1,3) - v(1,3,3) - x(3,3) = 2
MCLinV(3,1,3):  - x(1,3) - v(3,1,3) + x(3,3) = 9
MCLin1(1,3,3):  u(1,3,3) >= -0
MCLin2(1,3,3):  -99 y(3,1,3) - u(1,3,3) <= -0
MCLin3(1,3,3):  102 y(3,1,3) + v(1,3,3) - u(1,3,3) <= 102
MCLin4(1,3,3):  102 y(3,1,3) - u(1,3,3) >= -0
MCLin5(1,3,3):  -99 y(3,1,3) + v(1,3,3) - u(1,3,3) >= -99
MCLin6(3,1,3):  u(3,1,3) >= -0
MCLin7(3,1,3):  -113 y(1,3,3) - u(3,1,3) <= -0
MCLin8(3,1,3):  88 y(1,3,3) + v(3,1,3) - u(3,1,3) <= 88
MCLin9(3,1,3):  88 y(1,3,3) - u(3,1,3) >= -0
MCLin10(3,1,3):  -113 y(1,3,3) + v(3,1,3) - u(3,1,3) >= -113
MCLinV(1,4,0):  x(1,0) - v(1,4,0) - x(4,0) = 4
MCLinV(4,1,0):  - x(1,0) - v(4,1,0) + x(4,0) = 9
MCLin1(1,4,0):  u(1,4,0) >= -0
MCLin2(1,4,0):  -108 y(4,1,0) - u(1,4,0) <= -0
MCLin3(1,4,0):  91 y(4,1,0) + v(1,4,0) - u(1,4,0) <= 91
MCLin4(1,4,0):  91 y(4,1,0) - u(1,4,0) >= -0
MCLin5(1,4,0):  -108 y(4,1,0) + v(1,4,0) - u(1,4,0) >= -108
MCLin6(4,1,0):  u(4,1,0) >= -0
MCLin7(4,1,0):  -104 y(1,4,0) - u(4,1,0) <= -0
MCLin8(4,1,0):  95 y(1,4,0) + v(4,1,0) - u(4,1,0) <= 95
MCLin9(4,1,0):  95 y(1,4,0) - u(4,1,0) >= -0
MCLin10(4,1,0):  -104 y(1,4,0) + v(4,1,0) - u(4,1,0) >= -104
MCLinV(1,4,1):  x(1,1) - v(1,4,1) - x(4,1) = 10
MCLinV(4,1,1):  - x(1,1) - v(4,1,1) + x(4,1) = 2
MCLin1(1,4,1):  u(1,4,1) >= -0
MCLin2(1,4,1):  -87 y(4,1,1) - u(1,4,1) <= -0
MCLin3(1,4,1):  112 y(4,1,1) + v(1,4,1) - u(1,4,1) <= 112
MCLin4(1,4,1):  112 y(4,1,1) - u(1,4,1) >= -0
MCLin5(1,4,1):  -87 y(4,1,1) + v(1,4,1) - u(1,4,1) >= -87
MCLin6(4,1,1):  u(4,1,1) >= -0
MCLin7(4,1,1):  -124 y(1,4,1) - u(4,1,1) <= -0
MCLin8(4,1,1):  75 y(1,4,1) + v(4,1,1) - u(4,1,1) <= 75
MCLin9(4,1,1):  75 y(1,4,1) - u(4,1,1) >= -0
MCLin10(4,1,1):  -124 y(1,4,1) + v(4,1,1) - u(4,1,1) >= -124
MCLinV(1,4,2):  x(1,2) - v(1,4,2) - x(4,2) = 2
MCLinV(4,1,2):  - x(1,2) - v(4,1,2) + x(4,2) = 7
MCLin1(1,4,2):  u(1,4,2) >= -0
MCLin2(1,4,2):  -114 y(4,1,2) - u(1,4,2) <= -0
MCLin3(1,4,2):  85 y(4,1,2) + v(1,4,2) - u(1,4,2) <= 85
MCLin4(1,4,2):  85 y(4,1,2) - u(1,4,2) >= -0
MCLin5(1,4,2):  -114 y(4,1,2) + v(1,4,2) - u(1,4,2) >= -114
MCLin6(4,1,2):  u(4,1,2) >= -0
MCLin7(4,1,2):  -94 y(1,4,2) - u(4,1,2) <= -0
MCLin8(4,1,2):  105 y(1,4,2) + v(4,1,2) - u(4,1,2) <= 105
MCLin9(4,1,2):  105 y(1,4,2) - u(4,1,2) >= -0
MCLin10(4,1,2):  -94 y(1,4,2) + v(4,1,2) - u(4,1,2) >= -94
MCLinV(1,4,3):  x(1,3) - v(1,4,3) - x(4,3) = 6
MCLinV(4,1,3):  - x(1,3) - v(4,1,3) + x(4,3) = 9
MCLin1(1,4,3):  u(1,4,3) >= -0
MCLin2(1,4,3):  -113 y(4,1,3) - u(1,4,3) <= -0
MCLin3(1,4,3):  86 y(4,1,3) + v(1,4,3) - u(1,4,3) <= 86
MCLin4(1,4,3):  86 y(4,1,3) - u(1,4,3) >= -0
MCLin5(1,4,3):  -113 y(4,1,3) + v(1,4,3) - u(1,4,3) >= -113
MCLin6(4,1,3):  u(4,1,3) >= -0
MCLin7(4,1,3):  -101 y(1,4,3) - u(4,1,3) <= -0
MCLin8(4,1,3):  98 y(1,4,3) + v(4,1,3) - u(4,1,3) <= 98
MCLin9(4,1,3):  98 y(1,4,3) - u(4,1,3) >= -0
MCLin10(4,1,3):  -101 y(1,4,3) + v(4,1,3) - u(4,1,3) >= -101
MCLinV(2,3,0):  x(2,0) - v(2,3,0) - x(3,0) = 7
MCLinV(3,2,0):  - x(2,0) - v(3,2,0) + x(3,0) = 10
MCLin1(2,3,0):  u(2,3,0) >= -0
MCLin2(2,3,0):  -111 y(3,2,0) - u(2,3,0) <= -0
MCLin3(2,3,0):  94 y(3,2,0) + v(2,3,0) - u(2,3,0) <= 94
MCLin4(2,3,0):  94 y(3,2,0) - u(2,3,0) >= -0
MCLin5(2,3,0):  -111 y(3,2,0) + v(2,3,0) - u(2,3,0) >= -111
MCLin6(3,2,0):  u(3,2,0) >= -0
MCLin7(3,2,0):  -111 y(2,3,0) - u(3,2,0) <= -0
MCLin8(3,2,0):  94 y(2,3,0) + v(3,2,0) - u(3,2,0) <= 94
MCLin9(3,2,0):  94 y(2,3,0) - u(3,2,0) >= -0
MCLin10(3,2,0):  -111 y(2,3,0) + v(3,2,0) - u(3,2,0) >= -111
MCLinV(2,3,1):  x(2,1) - v(2,3,1) - x(3,1) = 10
MCLinV(3,2,1):  - x(2,1) - v(3,2,1) + x(3,1) = 2
MCLin1(2,3,1):  u(2,3,1) >= -0
MCLin2(2,3,1):  -106 y(3,2,1) - u(2,3,1) <= -0
MCLin3(2,3,1):  99 y(3,2,1) + v(2,3,1) - u(2,3,1) <= 99
MCLin4(2,3,1):  99 y(3,2,1) - u(2,3,1) >= -0
MCLin5(2,3,1):  -106 y(3,2,1) + v(2,3,1) - u(2,3,1) >= -106
MCLin6(3,2,1):  u(3,2,1) >= -0
MCLin7(3,2,1):  -111 y(2,3,1) - u(3,2,1) <= -0
MCLin8(3,2,1):  94 y(2,3,1) + v(3,2,1) - u(3,2,1) <= 94
MCLin9(3,2,1):  94 y(2,3,1) - u(3,2,1) >= -0
MCLin10(3,2,1):  -111 y(2,3,1) + v(3,2,1) - u(3,2,1) >= -111
MCLinV(2,3,2):  x(2,2) - v(2,3,2) - x(3,2) = 1
MCLinV(3,2,2):  - x(2,2) - v(3,2,2) + x(3,2) = 1
MCLin1(2,3,2):  u(2,3,2) >= -0
MCLin2(2,3,2):  -102 y(3,2,2) - u(2,3,2) <= -0
MCLin3(2,3,2):  103 y(3,2,2) + v(2,3,2) - u(2,3,2) <= 103
MCLin4(2,3,2):  103 y(3,2,2) - u(2,3,2) >= -0
MCLin5(2,3,2):  -102 y(3,2,2) + v(2,3,2) - u(2,3,2) >= -102
MCLin6(3,2,2):  u(3,2,2) >= -0
MCLin7(3,2,2):  -105 y(2,3,2) - u(3,2,2) <= -0
MCLin8(3,2,2):  100 y(2,3,2) + v(3,2,2) - u(3,2,2) <= 100
MCLin9(3,2,2):  100 y(2,3,2) - u(3,2,2) >= -0
MCLin10(3,2,2):  -105 y(2,3,2) + v(3,2,2) - u(3,2,2) >= -105
MCLinV(2,3,3):  x(2,3) - v(2,3,3) - x(3,3) = 2
MCLinV(3,2,3):  - x(2,3) - v(3,2,3) + x(3,3) = 10
MCLin1(2,3,3):  u(2,3,3) >= -0
MCLin2(2,3,3):  -106 y(3,2,3) - u(2,3,3) <= -0
MCLin3(2,3,3):  99 y(3,2,3) + v(2,3,3) - u(2,3,3) <= 99
MCLin4(2,3,3):  99 y(3,2,3) - u(2,3,3) >= -0
MCLin5(2,3,3):  -106 y(3,2,3) + v(2,3,3) - u(2,3,3) >= -106
MCLin6(3,2,3):  u(3,2,3) >= -0
MCLin7(3,2,3):  -111 y(2,3,3) - u(3,2,3) <= -0
MCLin8(3,2,3):  94 y(2,3,3) + v(3,2,3) - u(3,2,3) <= 94
MCLin9(3,2,3):  94 y(2,3,3) - u(3,2,3) >= -0
MCLin10(3,2,3):  -111 y(2,3,3) + v(3,2,3) - u(3,2,3) >= -111
MCLinV(2,4,0):  x(2,0) - v(2,4,0) - x(4,0) = 4
MCLinV(4,2,0):  - x(2,0) - v(4,2,0) + x(4,0) = 10
MCLin1(2,4,0):  u(2,4,0) >= -0
MCLin2(2,4,0):  -112 y(4,2,0) - u(2,4,0) <= -0
MCLin3(2,4,0):  91 y(4,2,0) + v(2,4,0) - u(2,4,0) <= 91
MCLin4(2,4,0):  91 y(4,2,0) - u(2,4,0) >= -0
MCLin5(2,4,0):  -112 y(4,2,0) + v(2,4,0) - u(2,4,0) >= -112
MCLin6(4,2,0):  u(4,2,0) >= -0
MCLin7(4,2,0):  -105 y(2,4,0) - u(4,2,0) <= -0
MCLin8(4,2,0):  98 y(2,4,0) + v(4,2,0) - u(4,2,0) <= 98
MCLin9(4,2,0):  98 y(2,4,0) - u(4,2,0) >= -0
MCLin10(4,2,0):  -105 y(2,4,0) + v(4,2,0) - u(4,2,0) >= -105
MCLinV(2,4,1):  x(2,1) - v(2,4,1) - x(4,1) = 10
MCLinV(4,2,1):  - x(2,1) - v(4,2,1) + x(4,1) = 2
MCLin1(2,4,1):  u(2,4,1) >= -0
MCLin2(2,4,1):  -102 y(4,2,1) - u(2,4,1) <= -0
MCLin3(2,4,1):  101 y(4,2,1) + v(2,4,1) - u(2,4,1) <= 101
MCLin4(2,4,1):  101 y(4,2,1) - u(2,4,1) >= -0
MCLin5(2,4,1):  -102 y(4,2,1) + v(2,4,1) - u(2,4,1) >= -102
MCLin6(4,2,1):  u(4,2,1) >= -0
MCLin7(4,2,1):  -113 y(2,4,1) - u(4,2,1) <= -0
MCLin8(4,2,1):  90 y(2,4,1) + v(4,2,1) - u(4,2,1) <= 90
MCLin9(4,2,1):  90 y(2,4,1) - u(4,2,1) >= -0
MCLin10(4,2,1):  -113 y(2,4,1) + v(4,2,1) - u(4,2,1) >= -113
MCLinV(2,4,2):  x(2,2) - v(2,4,2) - x(4,2) = 2
MCLinV(4,2,2):  - x(2,2) - v(4,2,2) + x(4,2) = 1
MCLin1(2,4,2):  u(2,4,2) >= -0
MCLin2(2,4,2):  -92 y(4,2,2) - u(2,4,2) <= -0
MCLin3(2,4,2):  111 y(4,2,2) + v(2,4,2) - u(2,4,2) <= 111
MCLin4(2,4,2):  111 y(4,2,2) - u(2,4,2) >= -0
MCLin5(2,4,2):  -92 y(4,2,2) + v(2,4,2) - u(2,4,2) >= -92
MCLin6(4,2,2):  u(4,2,2) >= -0
MCLin7(4,2,2):  -114 y(2,4,2) - u(4,2,2) <= -0
MCLin8(4,2,2):  89 y(2,4,2) + v(4,2,2) - u(4,2,2) <= 89
MCLin9(4,2,2):  89 y(2,4,2) - u(4,2,2) >= -0
MCLin10(4,2,2):  -114 y(2,4,2) + v(4,2,2) - u(4,2,2) >= -114
MCLinV(2,4,3):  x(2,3) - v(2,4,3) - x(4,3) = 6
MCLinV(4,2,3):  - x(2,3) - v(4,2,3) + x(4,3) = 10
MCLin1(2,4,3):  u(2,4,3) >= -0
MCLin2(2,4,3):  -120 y(4,2,3) - u(2,4,3) <= -0
MCLin3(2,4,3):  83 y(4,2,3) + v(2,4,3) - u(2,4,3) <= 83
MCLin4(2,4,3):  83 y(4,2,3) - u(2,4,3) >= -0
MCLin5(2,4,3):  -120 y(4,2,3) + v(2,4,3) - u(2,4,3) >= -120
MCLin6(4,2,3):  u(4,2,3) >= -0
MCLin7(4,2,3):  -99 y(2,4,3) - u(4,2,3) <= -0
MCLin8(4,2,3):  104 y(2,4,3) + v(4,2,3) - u(4,2,3) <= 104
MCLin9(4,2,3):  104 y(2,4,3) - u(4,2,3) >= -0
MCLin10(4,2,3):  -99 y(2,4,3) + v(4,2,3) - u(4,2,3) >= -99
MCLinV(3,4,0):  x(3,0) - v(3,4,0) - x(4,0) = 4
MCLinV(4,3,0):  - x(3,0) - v(4,3,0) + x(4,0) = 7
MCLin1(3,4,0):  u(3,4,0) >= -0
MCLin2(3,4,0):  -112 y(4,3,0) - u(3,4,0) <= -0
MCLin3(3,4,0):  94 y(4,3,0) + v(3,4,0) - u(3,4,0) <= 94
MCLin4(3,4,0):  94 y(4,3,0) - u(3,4,0) >= -0
MCLin5(3,4,0):  -112 y(4,3,0) + v(3,4,0) - u(3,4,0) >= -112
MCLin6(4,3,0):  u(4,3,0) >= -0
MCLin7(4,3,0):  -105 y(3,4,0) - u(4,3,0) <= -0
MCLin8(4,3,0):  101 y(3,4,0) + v(4,3,0) - u(4,3,0) <= 101
MCLin9(4,3,0):  101 y(3,4,0) - u(4,3,0) >= -0
MCLin10(4,3,0):  -105 y(3,4,0) + v(4,3,0) - u(4,3,0) >= -105
MCLinV(3,4,1):  x(3,1) - v(3,4,1) - x(4,1) = 10
MCLinV(4,3,1):  - x(3,1) - v(4,3,1) + x(4,1) = 10
MCLin1(3,4,1):  u(3,4,1) >= -0
MCLin2(3,4,1):  -110 y(4,3,1) - u(3,4,1) <= -0
MCLin3(3,4,1):  96 y(4,3,1) + v(3,4,1) - u(3,4,1) <= 96
MCLin4(3,4,1):  96 y(4,3,1) - u(3,4,1) >= -0
MCLin5(3,4,1):  -110 y(4,3,1) + v(3,4,1) - u(3,4,1) >= -110
MCLin6(4,3,1):  u(4,3,1) >= -0
MCLin7(4,3,1):  -116 y(3,4,1) - u(4,3,1) <= -0
MCLin8(4,3,1):  90 y(3,4,1) + v(4,3,1) - u(4,3,1) <= 90
MCLin9(4,3,1):  90 y(3,4,1) - u(4,3,1) >= -0
MCLin10(4,3,1):  -116 y(3,4,1) + v(4,3,1) - u(4,3,1) >= -116
MCLinV(3,4,2):  x(3,2) - v(3,4,2) - x(4,2) = 2
MCLinV(4,3,2):  - x(3,2) - v(4,3,2) + x(4,2) = 1
MCLin1(3,4,2):  u(3,4,2) >= -0
MCLin2(3,4,2):  -95 y(4,3,2) - u(3,4,2) <= -0
MCLin3(3,4,2):  111 y(4,3,2) + v(3,4,2) - u(3,4,2) <= 111
MCLin4(3,4,2):  111 y(4,3,2) - u(3,4,2) >= -0
MCLin5(3,4,2):  -95 y(4,3,2) + v(3,4,2) - u(3,4,2) >= -95
MCLin6(4,3,2):  u(4,3,2) >= -0
MCLin7(4,3,2):  -114 y(3,4,2) - u(4,3,2) <= -0
MCLin8(4,3,2):  92 y(3,4,2) + v(4,3,2) - u(4,3,2) <= 92
MCLin9(4,3,2):  92 y(3,4,2) - u(4,3,2) >= -0
MCLin10(4,3,2):  -114 y(3,4,2) + v(4,3,2) - u(4,3,2) >= -114
MCLinV(3,4,3):  x(3,3) - v(3,4,3) - x(4,3) = 6
MCLinV(4,3,3):  - x(3,3) - v(4,3,3) + x(4,3) = 2
MCLin1(3,4,3):  u(3,4,3) >= -0
MCLin2(3,4,3):  -120 y(4,3,3) - u(3,4,3) <= -0
MCLin3(3,4,3):  86 y(4,3,3) + v(3,4,3) - u(3,4,3) <= 86
MCLin4(3,4,3):  86 y(4,3,3) - u(3,4,3) >= -0
MCLin5(3,4,3):  -120 y(4,3,3) + v(3,4,3) - u(3,4,3) >= -120
MCLin6(4,3,3):  u(4,3,3) >= -0
MCLin7(4,3,3):  -94 y(3,4,3) - u(4,3,3) <= -0
MCLin8(4,3,3):  112 y(3,4,3) + v(4,3,3) - u(4,3,3) <= 112
MCLin9(4,3,3):  112 y(3,4,3) - u(4,3,3) >= -0
MCLin10(4,3,3):  -94 y(3,4,3) + v(4,3,3) - u(4,3,3) >= -94
makespan(0):  C - x(0,0) >= 6
makespan(1):  C - x(1,1) >= 2
makespan(2):  C - x(2,2) >= 1
makespan(3):  C - x(3,2) >= 1
makespan(4):  C - x(4,0) >= 4
cut_clique_0(01234,0):  0.00639 x(0,0) + 0.00958 x(1,0) + 0.01065 x(2,0) + 0.00745 x(3,0) + 0.00426 x(4,0) >= 1
cut_clique_0(1234,0):  0.01310 x(1,0) + 0.01456 x(2,0) + 0.01019 x(3,0) + 0.00582 x(4,0) >= 1
cut_clique_0(0134,0):  0.01529 x(0,0) + 0.01019 x(1,0) + 0.00156 x(3,0) + 0.01355 x(4,0) >= 1
cut_clique_0(0234,0):  0.01473 x(0,0) + 0.01146 x(2,0) + 0.01048 x(3,0) + 0.00655 x(4,0) >= 1
cut_clique_0(0123,0):  0.00786 x(0,0) + 0.01180 x(1,0) + 0.01311 x(2,0) + 0.00917 x(3,0) >= 1
cut_clique_0(0124,0):  0.00826 x(0,0) + 0.01239 x(1,0) + 0.00564 x(2,0) + 0.01459 x(4,0) >= 1
cut_clique_0(234,0):  0.01279 x(2,0) + 0.00242 x(3,0) + 0.03562 x(4,0) >= 1
cut_clique_0(134,0):  0.02438 x(1,0) + 0.00242 x(3,0) + 0.02159 x(4,0) >= 1
cut_clique_0(123,0):  0.01682 x(1,0) + 0.01869 x(2,0) + 0.01308 x(3,0) >= 1
cut_clique_0(124,0):  0.01906 x(1,0) + 0.00793 x(2,0) + 0.01860 x(4,0) >= 1
cut_clique_0(013,0):  0.03086 x(0,0) + 0.00475 x(1,0) + 0.00380 x(3,0) >= 1
cut_clique_0(014,0):  0.01823 x(0,0) + 0.00911 x(1,0) + 0.01302 x(4,0) >= 1
cut_clique_0(034,0):  0.00696 x(0,0) + 0.00209 x(3,0) + 0.04178 x(4,0) >= 1
cut_clique_0(012,0):  0.02137 x(0,0) + 0.01068 x(1,0) + 0.00855 x(2,0) >= 1
cut_clique_0(023,0):  0.02766 x(0,0) + 0.00638 x(2,0) + 0.00638 x(3,0) >= 1
cut_clique_0(024,0):  0.01008 x(0,0) + 0.00605 x(2,0) + 0.03024 x(4,0) >= 1
cut_clique_0(01234,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(0234,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(0123,1):  0.04348 x(0,1) + 0.02174 x(2,1) + 0.06522 x(3,1) >= 1
cut_clique_0(1234,1):  0.01429 x(2,1) + 0.07143 x(3,1) + 0.04898 x(4,1) >= 1
cut_clique_0(0134,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(0124,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(034,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(023,1):  0.04348 x(0,1) + 0.02174 x(2,1) + 0.06522 x(3,1) >= 1
cut_clique_0(024,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(234,1):  0.01429 x(2,1) + 0.07143 x(3,1) + 0.04898 x(4,1) >= 1
cut_clique_0(013,1):  0.06481 x(0,1) + 0.11111 x(3,1) >= 1
cut_clique_0(014,1):  0.10000 x(0,1) + 0.11111 x(4,1) >= 1
cut_clique_0(134,1):  0.10000 x(3,1) + 0.06667 x(4,1) >= 1
cut_clique_0(01234,2):  0.02299 x(0,2) + 0.01149 x(2,2) + 0.01149 x(3,2) >= 1
cut_clique_0(0234,2):  0.02299 x(0,2) + 0.01149 x(2,2) + 0.01149 x(3,2) >= 1
cut_clique_0(0134,2):  0.03553 x(0,2) + 0.01523 x(3,2) >= 1
cut_clique_0(0123,2):  0.02299 x(0,2) + 0.01149 x(2,2) + 0.01149 x(3,2) >= 1
cut_clique_0(023,2):  0.02299 x(0,2) + 0.01149 x(2,2) + 0.01149 x(3,2) >= 1
cut_clique_0(034,2):  0.03553 x(0,2) + 0.01523 x(3,2) >= 1
cut_clique_0(013,2):  0.03553 x(0,2) + 0.01523 x(3,2) >= 1
cut_clique_0(01234,3):  0.01746 x(0,3) + 0.01964 x(1,3) + 0.00896 x(2,3) + 0.00359 x(3,3) + 0.01309 x(4,3) >= 1
cut_clique_0(0134,3):  0.02541 x(0,3) + 0.01945 x(1,3) + 0.01906 x(4,3) >= 1
cut_clique_0(0123,3):  0.02820 x(0,3) + 0.03172 x(1,3) + 0.01362 x(2,3) + 0.00545 x(3,3) >= 1
cut_clique_0(0124,3):  0.01980 x(0,3) + 0.02228 x(1,3) + 0.00569 x(2,3) + 0.01485 x(4,3) >= 1
cut_clique_0(0234,3):  0.03636 x(0,3) + 0.00909 x(2,3) + 0.00606 x(3,3) + 0.02727 x(4,3) >= 1
cut_clique_0(1234,3):  0.03267 x(1,3) + 0.01149 x(2,3) + 0.00230 x(3,3) + 0.03241 x(4,3) >= 1
cut_clique_0(014,3):  0.02541 x(0,3) + 0.01945 x(1,3) + 0.01906 x(4,3) >= 1
cut_clique_0(023,3):  0.08065 x(0,3) + 0.01613 x(2,3) + 0.01075 x(3,3) >= 1
cut_clique_0(013,3):  0.04785 x(0,3) + 0.03349 x(1,3) >= 1
cut_clique_0(034,3):  0.02525 x(0,3) + 0.04545 x(4,3) >= 1
cut_clique_0(134,3):  0.01961 x(1,3) + 0.05392 x(4,3) >= 1
cut_clique_0(012,3):  0.03306 x(0,3) + 0.03719 x(1,3) + 0.00878 x(2,3) >= 1
cut_clique_0(024,3):  0.02941 x(0,3) + 0.00298 x(2,3) + 0.03922 x(4,3) >= 1
cut_clique_0(124,3):  0.03153 x(1,3) + 0.00921 x(2,3) + 0.03604 x(4,3) >= 1
cut_clique_1(01234,1):  0.02663 x(0,1) + 0.00592 x(2,1) + 0.02959 x(3,1) + 0.02959 x(4,1) >= 1
cut_clique_1(0234,1):  0.02663 x(0,1) + 0.00592 x(2,1) + 0.02959 x(3,1) + 0.02959 x(4,1) >= 1
cut_clique_1(0134,1):  0.03214 x(0,1) + 0.03571 x(3,1) + 0.03571 x(4,1) >= 1
cut_clique_1(034,1):  0.03214 x(0,1) + 0.03571 x(3,1) + 0.03571 x(4,1) >= 1
cut_clique_1(01234,3):  0.02835 x(0,3) + 0.01653 x(1,3) + 0.00554 x(2,3) + 0.00367 x(3,3) + 0.01102 x(4,3) >= 1
cut_clique_0:  0.00357 x(0,0) + 0.00535 x(1,0) + 0.00594 x(2,0) + 0.00416 x(3,0) + 0.00238 x(4,0) + 0.01937 x(0,1) + 0.00998 x(3,1) + 0.01854 x(4,1) >= 1
Bounds
 0 <= C <= 124
 26 <= x(0,0) <= 118
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -96 <= v(0,1,0) <= 93
 -108 <= v(1,0,0) <= 81
 -1240 <= u(0,1,0)
 -1240 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -97 <= v(0,2,0) <= 96
 -112 <= v(2,0,0) <= 81
 -1240 <= u(0,2,0)
 -1240 <= u(2,0,0)
 0 <= y(0,3,0) <= 1
 0 <= y(3,0,0) <= 1
 -97 <= v(0,3,0) <= 99
 -112 <= v(3,0,0) <= 84
 -1240 <= u(0,3,0)
 -1240 <= u(3,0,0)
 0 <= y(0,4,0) <= 1
 0 <= y(4,0,0) <= 1
 -98 <= v(0,4,0) <= 96
 -106 <= v(4,0,0) <= 88
 -1240 <= u(0,4,0)
 -1240 <= u(4,0,0)
 16 <= x(1,0) <= 113
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -107 <= v(1,2,0) <= 91
 -110 <= v(2,1,0) <= 88
 -1240 <= u(1,2,0)
 -1240 <= u(2,1,0)
 0 <= y(1,3,0) <= 1
 0 <= y(3,1,0) <= 1
 -107 <= v(1,3,0) <= 94
 -110 <= v(3,1,0) <= 91
 -1240 <= u(1,3,0)
 -1240 <= u(3,1,0)
 0 <= y(1,4,0) <= 1
 0 <= y(4,1,0) <= 1
 -108 <= v(1,4,0) <= 91
 -104 <= v(4,1,0) <= 95
 -1240 <= u(1,4,0)
 -1240 <= u(4,1,0)
 12 <= x(2,0) <= 113
 0 <= y(2,3,0) <= 1
 0 <= y(3,2,0) <= 1
 -111 <= v(2,3,0) <= 94
 -111 <= v(3,2,0) <= 94
 -1240 <= u(2,3,0)
 -1240 <= u(3,2,0)
 0 <= y(2,4,0) <= 1
 0 <= y(4,2,0) <= 1
 -112 <= v(2,4,0) <= 91
 -105 <= v(4,2,0) <= 98
 -1240 <= u(2,4,0)
 -1240 <= u(4,2,0)
 12 <= x(3,0) <= 116
 0 <= y(3,4,0) <= 1
 0 <= y(4,3,0) <= 1
 -112 <= v(3,4,0) <= 94
 -105 <= v(4,3,0) <= 101
 -1240 <= u(3,4,0)
 -1240 <= u(4,3,0)
 18 <= x(4,0) <= 120
 0 <= x(0,1) <= 92
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -124 <= v(0,1,1) <= 65
 -76 <= v(1,0,1) <= 113
 -1240 <= u(0,1,1)
 -1240 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -113 <= v(0,2,1) <= 80
 -91 <= v(2,0,1) <= 102
 -1240 <= u(0,2,1)
 -1240 <= u(2,0,1)
 0 <= y(0,3,1) <= 1
 0 <= y(3,0,1) <= 1
 -116 <= v(0,3,1) <= 80
 -99 <= v(3,0,1) <= 97
 -1240 <= u(0,3,1)
 -1240 <= u(3,0,1)
 0 <= y(0,4,1) <= 1
 0 <= y(4,0,1) <= 1
 -112 <= v(0,4,1) <= 82
 -101 <= v(4,0,1) <= 93
 -1240 <= u(0,4,1)
 -1240 <= u(4,0,1)
 25 <= x(1,1) <= 122
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -88 <= v(1,2,1) <= 110
 -114 <= v(2,1,1) <= 84
 -1240 <= u(1,2,1)
 -1240 <= u(2,1,1)
 0 <= y(1,3,1) <= 1
 0 <= y(3,1,1) <= 1
 -91 <= v(1,3,1) <= 110
 -122 <= v(3,1,1) <= 79
 -1240 <= u(1,3,1)
 -1240 <= u(3,1,1)
 0 <= y(1,4,1) <= 1
 0 <= y(4,1,1) <= 1
 -87 <= v(1,4,1) <= 112
 -124 <= v(4,1,1) <= 75
 -1240 <= u(1,4,1)
 -1240 <= u(4,1,1)
 10 <= x(2,1) <= 111
 0 <= y(2,3,1) <= 1
 0 <= y(3,2,1) <= 1
 -106 <= v(2,3,1) <= 99
 -111 <= v(3,2,1) <= 94
 -1240 <= u(2,3,1)
 -1240 <= u(3,2,1)
 0 <= y(2,4,1) <= 1
 0 <= y(4,2,1) <= 1
 -102 <= v(2,4,1) <= 101
 -113 <= v(4,2,1) <= 90
 -1240 <= u(2,4,1)
 -1240 <= u(4,2,1)
 2 <= x(3,1) <= 106
 0 <= y(3,4,1) <= 1
 0 <= y(4,3,1) <= 1
 -110 <= v(3,4,1) <= 96
 -116 <= v(4,3,1) <= 90
 -1240 <= u(3,4,1)
 -1240 <= u(4,3,1)
 0 <= x(4,1) <= 102
 17 <= x(0,2) <= 109
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -87 <= v(0,1,2) <= 102
 -118 <= v(1,0,2) <= 71
 -1240 <= u(0,1,2)
 -1240 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -107 <= v(0,2,2) <= 86
 -96 <= v(2,0,2) <= 97
 -1240 <= u(0,2,2)
 -1240 <= u(2,0,2)
 0 <= y(0,3,2) <= 1
 0 <= y(3,0,2) <= 1
 -107 <= v(0,3,2) <= 89
 -99 <= v(3,0,2) <= 97
 -1240 <= u(0,3,2)
 -1240 <= u(3,0,2)
 0 <= y(0,4,2) <= 1
 0 <= y(4,0,2) <= 1
 -97 <= v(0,4,2) <= 97
 -108 <= v(4,0,2) <= 86
 -1240 <= u(0,4,2)
 -1240 <= u(4,0,2)
 0 <= x(1,2) <= 97
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -124 <= v(1,2,2) <= 74
 -82 <= v(2,1,2) <= 116
 -1240 <= u(1,2,2)
 -1240 <= u(2,1,2)
 0 <= y(1,3,2) <= 1
 0 <= y(3,1,2) <= 1
 -124 <= v(1,3,2) <= 77
 -85 <= v(3,1,2) <= 116
 -1240 <= u(1,3,2)
 -1240 <= u(3,1,2)
 0 <= y(1,4,2) <= 1
 0 <= y(4,1,2) <= 1
 -114 <= v(1,4,2) <= 85
 -94 <= v(4,1,2) <= 105
 -1240 <= u(1,4,2)
 -1240 <= u(4,1,2)
 22 <= x(2,2) <= 123
 0 <= y(2,3,2) <= 1
 0 <= y(3,2,2) <= 1
 -102 <= v(2,3,2) <= 103
 -105 <= v(3,2,2) <= 100
 -1240 <= u(2,3,2)
 -1240 <= u(3,2,2)
 0 <= y(2,4,2) <= 1
 0 <= y(4,2,2) <= 1
 -92 <= v(2,4,2) <= 111
 -114 <= v(4,2,2) <= 89
 -1240 <= u(2,4,2)
 -1240 <= u(4,2,2)
 19 <= x(3,2) <= 123
 0 <= y(3,4,2) <= 1
 0 <= y(4,3,2) <= 1
 -95 <= v(3,4,2) <= 111
 -114 <= v(4,3,2) <= 92
 -1240 <= u(3,4,2)
 -1240 <= u(4,3,2)
 10 <= x(4,2) <= 112
 9 <= x(0,3) <= 101
 0 <= y(0,1,3) <= 1
 0 <= y(1,0,3) <= 1
 -104 <= v(0,1,3) <= 85
 -102 <= v(1,0,3) <= 87
 -1240 <= u(0,1,3)
 -1240 <= u(1,0,3)
 0 <= y(0,2,3) <= 1
 0 <= y(2,0,3) <= 1
 -102 <= v(0,2,3) <= 91
 -109 <= v(2,0,3) <= 84
 -1240 <= u(0,2,3)
 -1240 <= u(2,0,3)
 0 <= y(0,3,3) <= 1
 0 <= y(3,0,3) <= 1
 -97 <= v(0,3,3) <= 99
 -109 <= v(3,0,3) <= 87
 -1240 <= u(0,3,3)
 -1240 <= u(3,0,3)
 0 <= y(0,4,3) <= 1
 0 <= y(4,0,3) <= 1
 -111 <= v(0,4,3) <= 83
 -97 <= v(4,0,3) <= 97
 -1240 <= u(0,4,3)
 -1240 <= u(4,0,3)
 7 <= x(1,3) <= 104
 0 <= y(1,2,3) <= 1
 0 <= y(2,1,3) <= 1
 -104 <= v(1,2,3) <= 94
 -113 <= v(2,1,3) <= 85
 -1240 <= u(1,2,3)
 -1240 <= u(2,1,3)
 0 <= y(1,3,3) <= 1
 0 <= y(3,1,3) <= 1
 -99 <= v(1,3,3) <= 102
 -113 <= v(3,1,3) <= 88
 -1240 <= u(1,3,3)
 -1240 <= u(3,1,3)
 0 <= y(1,4,3) <= 1
 0 <= y(4,1,3) <= 1
 -113 <= v(1,4,3) <= 86
 -101 <= v(4,1,3) <= 98
 -1240 <= u(1,4,3)
 -1240 <= u(4,1,3)
 0 <= x(2,3) <= 101
 0 <= y(2,3,3) <= 1
 0 <= y(3,2,3) <= 1
 -106 <= v(2,3,3) <= 99
 -111 <= v(3,2,3) <= 94
 -1240 <= u(2,3,3)
 -1240 <= u(3,2,3)
 0 <= y(2,4,3) <= 1
 0 <= y(4,2,3) <= 1
 -120 <= v(2,4,3) <= 83
 -99 <= v(4,2,3) <= 104
 -1240 <= u(2,4,3)
 -1240 <= u(4,2,3)
 0 <= x(3,3) <= 104
 0 <= y(3,4,3) <= 1
 0 <= y(4,3,3) <= 1
 -120 <= v(3,4,3) <= 86
 -94 <= v(4,3,3) <= 112
 -1240 <= u(3,4,3)
 -1240 <= u(4,3,3)
 12 <= x(4,3) <= 114
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
