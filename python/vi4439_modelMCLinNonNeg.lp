\ Model compact
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C
Subject To
 ord(0,1): - x(0,0) + x(0,3) >= 10
 ord(0,2): x(0,2) - x(0,3) >= 1
 ord(0,3): x(0,1) - x(0,2) >= 6
 ord(1,1): - x(1,0) + x(1,3) >= 1
 ord(1,2): x(1,1) - x(1,3) >= 7
 ord(1,3): - x(1,1) + x(1,2) >= 5
 ord(2,1): x(2,2) - x(2,3) >= 9
 ord(2,2): x(2,1) - x(2,2) >= 3
 ord(2,3): x(2,0) - x(2,1) >= 7
 ord(3,1): - x(3,1) + x(3,3) >= 7
 ord(3,2): x(3,0) - x(3,3) >= 3
 ord(3,3): - x(3,0) + x(3,2) >= 9
 triangle2(0,1,0): y(0,1,0) + y(1,0,0) = 1
 triangle2(0,1,1): y(0,1,1) + y(1,0,1) = 1
 triangle2(0,1,2): y(0,1,2) + y(1,0,2) = 1
 triangle2(0,1,3): y(0,1,3) + y(1,0,3) = 1
 triangle2(0,2,0): y(0,2,0) + y(2,0,0) = 1
 triangle2(0,2,1): y(0,2,1) + y(2,0,1) = 1
 triangle2(0,2,2): y(0,2,2) + y(2,0,2) = 1
 triangle2(0,2,3): y(0,2,3) + y(2,0,3) = 1
 triangle2(0,3,0): y(0,3,0) + y(3,0,0) = 1
 triangle2(0,3,1): y(0,3,1) + y(3,0,1) = 1
 triangle2(0,3,2): y(0,3,2) + y(3,0,2) = 1
 triangle2(0,3,3): y(0,3,3) + y(3,0,3) = 1
 triangle2(1,2,0): y(1,2,0) + y(2,1,0) = 1
 triangle2(1,2,1): y(1,2,1) + y(2,1,1) = 1
 triangle2(1,2,2): y(1,2,2) + y(2,1,2) = 1
 triangle2(1,2,3): y(1,2,3) + y(2,1,3) = 1
 triangle2(1,3,0): y(1,3,0) + y(3,1,0) = 1
 triangle2(1,3,1): y(1,3,1) + y(3,1,1) = 1
 triangle2(1,3,2): y(1,3,2) + y(3,1,2) = 1
 triangle2(1,3,3): y(1,3,3) + y(3,1,3) = 1
 triangle2(2,3,0): y(2,3,0) + y(3,2,0) = 1
 triangle2(2,3,1): y(2,3,1) + y(3,2,1) = 1
 triangle2(2,3,2): y(2,3,2) + y(3,2,2) = 1
 triangle2(2,3,3): y(2,3,3) + y(3,2,3) = 1
 MCLinV(0,1,0): - x(0,0) + v(0,1,0)p - v(0,1,0)n + x(1,0) = -1
 MCLinV(1,0,0): x(0,0) + v(1,0,0)p - v(1,0,0)n - x(1,0) = -10
 MCLin1(0,1,0): u(0,1,0)p - u(0,1,0)n >= -0
 MCLin2(0,1,0): 33 y(1,0,0) + u(0,1,0)p - u(0,1,0)n >= -0
 MCLin3(0,1,0): - 28 y(1,0,0) - v(0,1,0)p + v(0,1,0)n + u(0,1,0)p
   - u(0,1,0)n >= -28
 MCLin4(0,1,0): - 28 y(1,0,0) + u(0,1,0)p - u(0,1,0)n <= -0
 MCLin5(0,1,0): 33 y(1,0,0) - v(0,1,0)p + v(0,1,0)n + u(0,1,0)p - u(0,1,0)n
   <= 33
 MCLin6(1,0,0): u(1,0,0)p - u(1,0,0)n >= -0
 MCLin7(1,0,0): 39 y(0,1,0) + u(1,0,0)p - u(1,0,0)n >= -0
 MCLin8(1,0,0): - 22 y(0,1,0) - v(1,0,0)p + v(1,0,0)n + u(1,0,0)p
   - u(1,0,0)n >= -22
 MCLin9(1,0,0): - 22 y(0,1,0) + u(1,0,0)p - u(1,0,0)n <= -0
 MCLin10(1,0,0): 39 y(0,1,0) - v(1,0,0)p + v(1,0,0)n + u(1,0,0)p
   - u(1,0,0)n <= 39
 MCLinV(0,1,1): - x(0,1) + v(0,1,1)p - v(0,1,1)n + x(1,1) = -5
 MCLinV(1,0,1): x(0,1) + v(1,0,1)p - v(1,0,1)n - x(1,1) = -7
 MCLin1(0,1,1): u(0,1,1)p - u(0,1,1)n >= -0
 MCLin2(0,1,1): 28 y(1,0,1) + u(0,1,1)p - u(0,1,1)n >= -0
 MCLin3(0,1,1): - 33 y(1,0,1) - v(0,1,1)p + v(0,1,1)n + u(0,1,1)p
   - u(0,1,1)n >= -33
 MCLin4(0,1,1): - 33 y(1,0,1) + u(0,1,1)p - u(0,1,1)n <= -0
 MCLin5(0,1,1): 28 y(1,0,1) - v(0,1,1)p + v(0,1,1)n + u(0,1,1)p - u(0,1,1)n
   <= 28
 MCLin6(1,0,1): u(1,0,1)p - u(1,0,1)n >= -0
 MCLin7(1,0,1): 45 y(0,1,1) + u(1,0,1)p - u(1,0,1)n >= -0
 MCLin8(1,0,1): - 16 y(0,1,1) - v(1,0,1)p + v(1,0,1)n + u(1,0,1)p
   - u(1,0,1)n >= -16
 MCLin9(1,0,1): - 16 y(0,1,1) + u(1,0,1)p - u(1,0,1)n <= -0
 MCLin10(1,0,1): 45 y(0,1,1) - v(1,0,1)p + v(1,0,1)n + u(1,0,1)p
   - u(1,0,1)n <= 45
 MCLinV(0,1,2): - x(0,2) + v(0,1,2)p - v(0,1,2)n + x(1,2) = -8
 MCLinV(1,0,2): x(0,2) + v(1,0,2)p - v(1,0,2)n - x(1,2) = -6
 MCLin1(0,1,2): u(0,1,2)p - u(0,1,2)n >= -0
 MCLin2(0,1,2): 42 y(1,0,2) + u(0,1,2)p - u(0,1,2)n >= -0
 MCLin3(0,1,2): - 19 y(1,0,2) - v(0,1,2)p + v(0,1,2)n + u(0,1,2)p
   - u(0,1,2)n >= -19
 MCLin4(0,1,2): - 19 y(1,0,2) + u(0,1,2)p - u(0,1,2)n <= -0
 MCLin5(0,1,2): 42 y(1,0,2) - v(0,1,2)p + v(0,1,2)n + u(0,1,2)p - u(0,1,2)n
   <= 42
 MCLin6(1,0,2): u(1,0,2)p - u(1,0,2)n >= -0
 MCLin7(1,0,2): 33 y(0,1,2) + u(1,0,2)p - u(1,0,2)n >= -0
 MCLin8(1,0,2): - 28 y(0,1,2) - v(1,0,2)p + v(1,0,2)n + u(1,0,2)p
   - u(1,0,2)n >= -28
 MCLin9(1,0,2): - 28 y(0,1,2) + u(1,0,2)p - u(1,0,2)n <= -0
 MCLin10(1,0,2): 33 y(0,1,2) - v(1,0,2)p + v(1,0,2)n + u(1,0,2)p
   - u(1,0,2)n <= 33
 MCLinV(0,1,3): - x(0,3) + v(0,1,3)p - v(0,1,3)n + x(1,3) = -7
 MCLinV(1,0,3): x(0,3) + v(1,0,3)p - v(1,0,3)n - x(1,3) = -1
 MCLin1(0,1,3): u(0,1,3)p - u(0,1,3)n >= -0
 MCLin2(0,1,3): 30 y(1,0,3) + u(0,1,3)p - u(0,1,3)n >= -0
 MCLin3(0,1,3): - 31 y(1,0,3) - v(0,1,3)p + v(0,1,3)n + u(0,1,3)p
   - u(0,1,3)n >= -31
 MCLin4(0,1,3): - 31 y(1,0,3) + u(0,1,3)p - u(0,1,3)n <= -0
 MCLin5(0,1,3): 30 y(1,0,3) - v(0,1,3)p + v(0,1,3)n + u(0,1,3)p - u(0,1,3)n
   <= 30
 MCLin6(1,0,3): u(1,0,3)p - u(1,0,3)n >= -0
 MCLin7(1,0,3): 39 y(0,1,3) + u(1,0,3)p - u(1,0,3)n >= -0
 MCLin8(1,0,3): - 22 y(0,1,3) - v(1,0,3)p + v(1,0,3)n + u(1,0,3)p
   - u(1,0,3)n >= -22
 MCLin9(1,0,3): - 22 y(0,1,3) + u(1,0,3)p - u(1,0,3)n <= -0
 MCLin10(1,0,3): 39 y(0,1,3) - v(1,0,3)p + v(1,0,3)n + u(1,0,3)p
   - u(1,0,3)n <= 39
 MCLinV(0,2,0): - x(0,0) + v(0,2,0)p - v(0,2,0)n + x(2,0) = -1
 MCLinV(2,0,0): x(0,0) + v(2,0,0)p - v(2,0,0)n - x(2,0) = -10
 MCLin1(0,2,0): u(0,2,0)p - u(0,2,0)n >= -0
 MCLin2(0,2,0): 53 y(2,0,0) + u(0,2,0)p - u(0,2,0)n >= -0
 MCLin3(0,2,0): - 9 y(2,0,0) - v(0,2,0)p + v(0,2,0)n + u(0,2,0)p
   - u(0,2,0)n >= -9
 MCLin4(0,2,0): - 9 y(2,0,0) + u(0,2,0)p - u(0,2,0)n <= -0
 MCLin5(0,2,0): 53 y(2,0,0) - v(0,2,0)p + v(0,2,0)n + u(0,2,0)p - u(0,2,0)n
   <= 53
 MCLin6(2,0,0): u(2,0,0)p - u(2,0,0)n >= -0
 MCLin7(2,0,0): 20 y(0,2,0) + u(2,0,0)p - u(2,0,0)n >= -0
 MCLin8(2,0,0): - 42 y(0,2,0) - v(2,0,0)p + v(2,0,0)n + u(2,0,0)p
   - u(2,0,0)n >= -42
 MCLin9(2,0,0): - 42 y(0,2,0) + u(2,0,0)p - u(2,0,0)n <= -0
 MCLin10(2,0,0): 20 y(0,2,0) - v(2,0,0)p + v(2,0,0)n + u(2,0,0)p
   - u(2,0,0)n <= 20
 MCLinV(0,2,1): - x(0,1) + v(0,2,1)p - v(0,2,1)n + x(2,1) = -7
 MCLinV(2,0,1): x(0,1) + v(2,0,1)p - v(2,0,1)n - x(2,1) = -7
 MCLin1(0,2,1): u(0,2,1)p - u(0,2,1)n >= -0
 MCLin2(0,2,1): 35 y(2,0,1) + u(0,2,1)p - u(0,2,1)n >= -0
 MCLin3(0,2,1): - 27 y(2,0,1) - v(0,2,1)p + v(0,2,1)n + u(0,2,1)p
   - u(0,2,1)n >= -27
 MCLin4(0,2,1): - 27 y(2,0,1) + u(0,2,1)p - u(0,2,1)n <= -0
 MCLin5(0,2,1): 35 y(2,0,1) - v(0,2,1)p + v(0,2,1)n + u(0,2,1)p - u(0,2,1)n
   <= 35
 MCLin6(2,0,1): u(2,0,1)p - u(2,0,1)n >= -0
 MCLin7(2,0,1): 41 y(0,2,1) + u(2,0,1)p - u(2,0,1)n >= -0
 MCLin8(2,0,1): - 21 y(0,2,1) - v(2,0,1)p + v(2,0,1)n + u(2,0,1)p
   - u(2,0,1)n >= -21
 MCLin9(2,0,1): - 21 y(0,2,1) + u(2,0,1)p - u(2,0,1)n <= -0
 MCLin10(2,0,1): 41 y(0,2,1) - v(2,0,1)p + v(2,0,1)n + u(2,0,1)p
   - u(2,0,1)n <= 41
 MCLinV(0,2,2): - x(0,2) + v(0,2,2)p - v(0,2,2)n + x(2,2) = -3
 MCLinV(2,0,2): x(0,2) + v(2,0,2)p - v(2,0,2)n - x(2,2) = -6
 MCLin1(0,2,2): u(0,2,2)p - u(0,2,2)n >= -0
 MCLin2(0,2,2): 34 y(2,0,2) + u(0,2,2)p - u(0,2,2)n >= -0
 MCLin3(0,2,2): - 28 y(2,0,2) - v(0,2,2)p + v(0,2,2)n + u(0,2,2)p
   - u(0,2,2)n >= -28
 MCLin4(0,2,2): - 28 y(2,0,2) + u(0,2,2)p - u(0,2,2)n <= -0
 MCLin5(0,2,2): 34 y(2,0,2) - v(0,2,2)p + v(0,2,2)n + u(0,2,2)p - u(0,2,2)n
   <= 34
 MCLin6(2,0,2): u(2,0,2)p - u(2,0,2)n >= -0
 MCLin7(2,0,2): 37 y(0,2,2) + u(2,0,2)p - u(2,0,2)n >= -0
 MCLin8(2,0,2): - 25 y(0,2,2) - v(2,0,2)p + v(2,0,2)n + u(2,0,2)p
   - u(2,0,2)n >= -25
 MCLin9(2,0,2): - 25 y(0,2,2) + u(2,0,2)p - u(2,0,2)n <= -0
 MCLin10(2,0,2): 37 y(0,2,2) - v(2,0,2)p + v(2,0,2)n + u(2,0,2)p
   - u(2,0,2)n <= 37
 MCLinV(0,2,3): - x(0,3) + v(0,2,3)p - v(0,2,3)n + x(2,3) = -9
 MCLinV(2,0,3): x(0,3) + v(2,0,3)p - v(2,0,3)n - x(2,3) = -1
 MCLin1(0,2,3): u(0,2,3)p - u(0,2,3)n >= -0
 MCLin2(0,2,3): 32 y(2,0,3) + u(0,2,3)p - u(0,2,3)n >= -0
 MCLin3(0,2,3): - 30 y(2,0,3) - v(0,2,3)p + v(0,2,3)n + u(0,2,3)p
   - u(0,2,3)n >= -30
 MCLin4(0,2,3): - 30 y(2,0,3) + u(0,2,3)p - u(0,2,3)n <= -0
 MCLin5(0,2,3): 32 y(2,0,3) - v(0,2,3)p + v(0,2,3)n + u(0,2,3)p - u(0,2,3)n
   <= 32
 MCLin6(2,0,3): u(2,0,3)p - u(2,0,3)n >= -0
 MCLin7(2,0,3): 40 y(0,2,3) + u(2,0,3)p - u(2,0,3)n >= -0
 MCLin8(2,0,3): - 22 y(0,2,3) - v(2,0,3)p + v(2,0,3)n + u(2,0,3)p
   - u(2,0,3)n >= -22
 MCLin9(2,0,3): - 22 y(0,2,3) + u(2,0,3)p - u(2,0,3)n <= -0
 MCLin10(2,0,3): 40 y(0,2,3) - v(2,0,3)p + v(2,0,3)n + u(2,0,3)p
   - u(2,0,3)n <= 40
 MCLinV(0,3,0): - x(0,0) + v(0,3,0)p - v(0,3,0)n + x(3,0) = -9
 MCLinV(3,0,0): x(0,0) + v(3,0,0)p - v(3,0,0)n - x(3,0) = -10
 MCLin1(0,3,0): u(0,3,0)p - u(0,3,0)n >= -0
 MCLin2(0,3,0): 43 y(3,0,0) + u(0,3,0)p - u(0,3,0)n >= -0
 MCLin3(0,3,0): - 10 y(3,0,0) - v(0,3,0)p + v(0,3,0)n + u(0,3,0)p
   - u(0,3,0)n >= -10
 MCLin4(0,3,0): - 10 y(3,0,0) + u(0,3,0)p - u(0,3,0)n <= -0
 MCLin5(0,3,0): 43 y(3,0,0) - v(0,3,0)p + v(0,3,0)n + u(0,3,0)p - u(0,3,0)n
   <= 43
 MCLin6(3,0,0): u(3,0,0)p - u(3,0,0)n >= -0
 MCLin7(3,0,0): 29 y(0,3,0) + u(3,0,0)p - u(3,0,0)n >= -0
 MCLin8(3,0,0): - 24 y(0,3,0) - v(3,0,0)p + v(3,0,0)n + u(3,0,0)p
   - u(3,0,0)n >= -24
 MCLin9(3,0,0): - 24 y(0,3,0) + u(3,0,0)p - u(3,0,0)n <= -0
 MCLin10(3,0,0): 29 y(0,3,0) - v(3,0,0)p + v(3,0,0)n + u(3,0,0)p
   - u(3,0,0)n <= 29
 MCLinV(0,3,1): - x(0,1) + v(0,3,1)p - v(0,3,1)n + x(3,1) = -7
 MCLinV(3,0,1): x(0,1) + v(3,0,1)p - v(3,0,1)n - x(3,1) = -7
 MCLin1(0,3,1): u(0,3,1)p - u(0,3,1)n >= -0
 MCLin2(0,3,1): 14 y(3,0,1) + u(0,3,1)p - u(0,3,1)n >= -0
 MCLin3(0,3,1): - 39 y(3,0,1) - v(0,3,1)p + v(0,3,1)n + u(0,3,1)p
   - u(0,3,1)n >= -39
 MCLin4(0,3,1): - 39 y(3,0,1) + u(0,3,1)p - u(0,3,1)n <= -0
 MCLin5(0,3,1): 14 y(3,0,1) - v(0,3,1)p + v(0,3,1)n + u(0,3,1)p - u(0,3,1)n
   <= 14
 MCLin6(3,0,1): u(3,0,1)p - u(3,0,1)n >= -0
 MCLin7(3,0,1): 53 y(0,3,1) + u(3,0,1)p - u(3,0,1)n >= -0
 MCLin8(3,0,1): - v(3,0,1)p + v(3,0,1)n + u(3,0,1)p - u(3,0,1)n >= -0
 MCLin9(3,0,1): u(3,0,1)p - u(3,0,1)n <= -0
 MCLin10(3,0,1): 53 y(0,3,1) - v(3,0,1)p + v(3,0,1)n + u(3,0,1)p
   - u(3,0,1)n <= 53
 MCLinV(0,3,2): - x(0,2) + v(0,3,2)p - v(0,3,2)n + x(3,2) = -10
 MCLinV(3,0,2): x(0,2) + v(3,0,2)p - v(3,0,2)n - x(3,2) = -6
 MCLin1(0,3,2): u(0,3,2)p - u(0,3,2)n >= -0
 MCLin2(0,3,2): 42 y(3,0,2) + u(0,3,2)p - u(0,3,2)n >= -0
 MCLin3(0,3,2): - 11 y(3,0,2) - v(0,3,2)p + v(0,3,2)n + u(0,3,2)p
   - u(0,3,2)n >= -11
 MCLin4(0,3,2): - 11 y(3,0,2) + u(0,3,2)p - u(0,3,2)n <= -0
 MCLin5(0,3,2): 42 y(3,0,2) - v(0,3,2)p + v(0,3,2)n + u(0,3,2)p - u(0,3,2)n
   <= 42
 MCLin6(3,0,2): u(3,0,2)p - u(3,0,2)n >= -0
 MCLin7(3,0,2): 27 y(0,3,2) + u(3,0,2)p - u(3,0,2)n >= -0
 MCLin8(3,0,2): - 26 y(0,3,2) - v(3,0,2)p + v(3,0,2)n + u(3,0,2)p
   - u(3,0,2)n >= -26
 MCLin9(3,0,2): - 26 y(0,3,2) + u(3,0,2)p - u(3,0,2)n <= -0
 MCLin10(3,0,2): 27 y(0,3,2) - v(3,0,2)p + v(3,0,2)n + u(3,0,2)p
   - u(3,0,2)n <= 27
 MCLinV(0,3,3): - x(0,3) + v(0,3,3)p - v(0,3,3)n + x(3,3) = -3
 MCLinV(3,0,3): x(0,3) + v(3,0,3)p - v(3,0,3)n - x(3,3) = -1
 MCLin1(0,3,3): u(0,3,3)p - u(0,3,3)n >= -0
 MCLin2(0,3,3): 24 y(3,0,3) + u(0,3,3)p - u(0,3,3)n >= -0
 MCLin3(0,3,3): - 29 y(3,0,3) - v(0,3,3)p + v(0,3,3)n + u(0,3,3)p
   - u(0,3,3)n >= -29
 MCLin4(0,3,3): - 29 y(3,0,3) + u(0,3,3)p - u(0,3,3)n <= -0
 MCLin5(0,3,3): 24 y(3,0,3) - v(0,3,3)p + v(0,3,3)n + u(0,3,3)p - u(0,3,3)n
   <= 24
 MCLin6(3,0,3): u(3,0,3)p - u(3,0,3)n >= -0
 MCLin7(3,0,3): 33 y(0,3,3) + u(3,0,3)p - u(3,0,3)n >= -0
 MCLin8(3,0,3): - 20 y(0,3,3) - v(3,0,3)p + v(3,0,3)n + u(3,0,3)p
   - u(3,0,3)n >= -20
 MCLin9(3,0,3): - 20 y(0,3,3) + u(3,0,3)p - u(3,0,3)n <= -0
 MCLin10(3,0,3): 33 y(0,3,3) - v(3,0,3)p + v(3,0,3)n + u(3,0,3)p
   - u(3,0,3)n <= 33
 MCLinV(1,2,0): - x(1,0) + v(1,2,0)p - v(1,2,0)n + x(2,0) = -1
 MCLinV(2,1,0): x(1,0) + v(2,1,0)p - v(2,1,0)n - x(2,0) = -1
 MCLin1(1,2,0): u(1,2,0)p - u(1,2,0)n >= -0
 MCLin2(1,2,0): 53 y(2,1,0) + u(1,2,0)p - u(1,2,0)n >= -0
 MCLin3(1,2,0): - 12 y(2,1,0) - v(1,2,0)p + v(1,2,0)n + u(1,2,0)p
   - u(1,2,0)n >= -12
 MCLin4(1,2,0): - 12 y(2,1,0) + u(1,2,0)p - u(1,2,0)n <= -0
 MCLin5(1,2,0): 53 y(2,1,0) - v(1,2,0)p + v(1,2,0)n + u(1,2,0)p - u(1,2,0)n
   <= 53
 MCLin6(2,1,0): u(2,1,0)p - u(2,1,0)n >= -0
 MCLin7(2,1,0): 14 y(1,2,0) + u(2,1,0)p - u(2,1,0)n >= -0
 MCLin8(2,1,0): - 51 y(1,2,0) - v(2,1,0)p + v(2,1,0)n + u(2,1,0)p
   - u(2,1,0)n >= -51
 MCLin9(2,1,0): - 51 y(1,2,0) + u(2,1,0)p - u(2,1,0)n <= -0
 MCLin10(2,1,0): 14 y(1,2,0) - v(2,1,0)p + v(2,1,0)n + u(2,1,0)p
   - u(2,1,0)n <= 14
 MCLinV(1,2,1): - x(1,1) + v(1,2,1)p - v(1,2,1)n + x(2,1) = -7
 MCLinV(2,1,1): x(1,1) + v(2,1,1)p - v(2,1,1)n - x(2,1) = -5
 MCLin1(1,2,1): u(1,2,1)p - u(1,2,1)n >= -0
 MCLin2(1,2,1): 44 y(2,1,1) + u(1,2,1)p - u(1,2,1)n >= -0
 MCLin3(1,2,1): - 21 y(2,1,1) - v(1,2,1)p + v(1,2,1)n + u(1,2,1)p
   - u(1,2,1)n >= -21
 MCLin4(1,2,1): - 21 y(2,1,1) + u(1,2,1)p - u(1,2,1)n <= -0
 MCLin5(1,2,1): 44 y(2,1,1) - v(1,2,1)p + v(1,2,1)n + u(1,2,1)p - u(1,2,1)n
   <= 44
 MCLin6(2,1,1): u(2,1,1)p - u(2,1,1)n >= -0
 MCLin7(2,1,1): 33 y(1,2,1) + u(2,1,1)p - u(2,1,1)n >= -0
 MCLin8(2,1,1): - 32 y(1,2,1) - v(2,1,1)p + v(2,1,1)n + u(2,1,1)p
   - u(2,1,1)n >= -32
 MCLin9(2,1,1): - 32 y(1,2,1) + u(2,1,1)p - u(2,1,1)n <= -0
 MCLin10(2,1,1): 33 y(1,2,1) - v(2,1,1)p + v(2,1,1)n + u(2,1,1)p
   - u(2,1,1)n <= 33
 MCLinV(1,2,2): - x(1,2) + v(1,2,2)p - v(1,2,2)n + x(2,2) = -3
 MCLinV(2,1,2): x(1,2) + v(2,1,2)p - v(2,1,2)n - x(2,2) = -8
 MCLin1(1,2,2): u(1,2,2)p - u(1,2,2)n >= -0
 MCLin2(1,2,2): 32 y(2,1,2) + u(1,2,2)p - u(1,2,2)n >= -0
 MCLin3(1,2,2): - 33 y(2,1,2) - v(1,2,2)p + v(1,2,2)n + u(1,2,2)p
   - u(1,2,2)n >= -33
 MCLin4(1,2,2): - 33 y(2,1,2) + u(1,2,2)p - u(1,2,2)n <= -0
 MCLin5(1,2,2): 32 y(2,1,2) - v(1,2,2)p + v(1,2,2)n + u(1,2,2)p - u(1,2,2)n
   <= 32
 MCLin6(2,1,2): u(2,1,2)p - u(2,1,2)n >= -0
 MCLin7(2,1,2): 44 y(1,2,2) + u(2,1,2)p - u(2,1,2)n >= -0
 MCLin8(2,1,2): - 21 y(1,2,2) - v(2,1,2)p + v(2,1,2)n + u(2,1,2)p
   - u(2,1,2)n >= -21
 MCLin9(2,1,2): - 21 y(1,2,2) + u(2,1,2)p - u(2,1,2)n <= -0
 MCLin10(2,1,2): 44 y(1,2,2) - v(2,1,2)p + v(2,1,2)n + u(2,1,2)p
   - u(2,1,2)n <= 44
 MCLinV(1,2,3): - x(1,3) + v(1,2,3)p - v(1,2,3)n + x(2,3) = -9
 MCLinV(2,1,3): x(1,3) + v(2,1,3)p - v(2,1,3)n - x(2,3) = -7
 MCLin1(1,2,3): u(1,2,3)p - u(1,2,3)n >= -0
 MCLin2(1,2,3): 41 y(2,1,3) + u(1,2,3)p - u(1,2,3)n >= -0
 MCLin3(1,2,3): - 24 y(2,1,3) - v(1,2,3)p + v(1,2,3)n + u(1,2,3)p
   - u(1,2,3)n >= -24
 MCLin4(1,2,3): - 24 y(2,1,3) + u(1,2,3)p - u(1,2,3)n <= -0
 MCLin5(1,2,3): 41 y(2,1,3) - v(1,2,3)p + v(1,2,3)n + u(1,2,3)p - u(1,2,3)n
   <= 41
 MCLin6(2,1,3): u(2,1,3)p - u(2,1,3)n >= -0
 MCLin7(2,1,3): 40 y(1,2,3) + u(2,1,3)p - u(2,1,3)n >= -0
 MCLin8(2,1,3): - 25 y(1,2,3) - v(2,1,3)p + v(2,1,3)n + u(2,1,3)p
   - u(2,1,3)n >= -25
 MCLin9(2,1,3): - 25 y(1,2,3) + u(2,1,3)p - u(2,1,3)n <= -0
 MCLin10(2,1,3): 40 y(1,2,3) - v(2,1,3)p + v(2,1,3)n + u(2,1,3)p
   - u(2,1,3)n <= 40
 MCLinV(1,3,0): - x(1,0) + v(1,3,0)p - v(1,3,0)n + x(3,0) = -9
 MCLinV(3,1,0): x(1,0) + v(3,1,0)p - v(3,1,0)n - x(3,0) = -1
 MCLin1(1,3,0): u(1,3,0)p - u(1,3,0)n >= -0
 MCLin2(1,3,0): 43 y(3,1,0) + u(1,3,0)p - u(1,3,0)n >= -0
 MCLin3(1,3,0): - 13 y(3,1,0) - v(1,3,0)p + v(1,3,0)n + u(1,3,0)p
   - u(1,3,0)n >= -13
 MCLin4(1,3,0): - 13 y(3,1,0) + u(1,3,0)p - u(1,3,0)n <= -0
 MCLin5(1,3,0): 43 y(3,1,0) - v(1,3,0)p + v(1,3,0)n + u(1,3,0)p - u(1,3,0)n
   <= 43
 MCLin6(3,1,0): u(3,1,0)p - u(3,1,0)n >= -0
 MCLin7(3,1,0): 23 y(1,3,0) + u(3,1,0)p - u(3,1,0)n >= -0
 MCLin8(3,1,0): - 33 y(1,3,0) - v(3,1,0)p + v(3,1,0)n + u(3,1,0)p
   - u(3,1,0)n >= -33
 MCLin9(3,1,0): - 33 y(1,3,0) + u(3,1,0)p - u(3,1,0)n <= -0
 MCLin10(3,1,0): 23 y(1,3,0) - v(3,1,0)p + v(3,1,0)n + u(3,1,0)p
   - u(3,1,0)n <= 23
 MCLinV(1,3,1): - x(1,1) + v(1,3,1)p - v(1,3,1)n + x(3,1) = -7
 MCLinV(3,1,1): x(1,1) + v(3,1,1)p - v(3,1,1)n - x(3,1) = -5
 MCLin1(1,3,1): u(1,3,1)p - u(1,3,1)n >= -0
 MCLin2(1,3,1): 23 y(3,1,1) + u(1,3,1)p - u(1,3,1)n >= -0
 MCLin3(1,3,1): - 33 y(3,1,1) - v(1,3,1)p + v(1,3,1)n + u(1,3,1)p
   - u(1,3,1)n >= -33
 MCLin4(1,3,1): - 33 y(3,1,1) + u(1,3,1)p - u(1,3,1)n <= -0
 MCLin5(1,3,1): 23 y(3,1,1) - v(1,3,1)p + v(1,3,1)n + u(1,3,1)p - u(1,3,1)n
   <= 23
 MCLin6(3,1,1): u(3,1,1)p - u(3,1,1)n >= -0
 MCLin7(3,1,1): 45 y(1,3,1) + u(3,1,1)p - u(3,1,1)n >= -0
 MCLin8(3,1,1): - 11 y(1,3,1) - v(3,1,1)p + v(3,1,1)n + u(3,1,1)p
   - u(3,1,1)n >= -11
 MCLin9(3,1,1): - 11 y(1,3,1) + u(3,1,1)p - u(3,1,1)n <= -0
 MCLin10(3,1,1): 45 y(1,3,1) - v(3,1,1)p + v(3,1,1)n + u(3,1,1)p
   - u(3,1,1)n <= 45
 MCLinV(1,3,2): - x(1,2) + v(1,3,2)p - v(1,3,2)n + x(3,2) = -10
 MCLinV(3,1,2): x(1,2) + v(3,1,2)p - v(3,1,2)n - x(3,2) = -8
 MCLin1(1,3,2): u(1,3,2)p - u(1,3,2)n >= -0
 MCLin2(1,3,2): 40 y(3,1,2) + u(1,3,2)p - u(1,3,2)n >= -0
 MCLin3(1,3,2): - 16 y(3,1,2) - v(1,3,2)p + v(1,3,2)n + u(1,3,2)p
   - u(1,3,2)n >= -16
 MCLin4(1,3,2): - 16 y(3,1,2) + u(1,3,2)p - u(1,3,2)n <= -0
 MCLin5(1,3,2): 40 y(3,1,2) - v(1,3,2)p + v(1,3,2)n + u(1,3,2)p - u(1,3,2)n
   <= 40
 MCLin6(3,1,2): u(3,1,2)p - u(3,1,2)n >= -0
 MCLin7(3,1,2): 34 y(1,3,2) + u(3,1,2)p - u(3,1,2)n >= -0
 MCLin8(3,1,2): - 22 y(1,3,2) - v(3,1,2)p + v(3,1,2)n + u(3,1,2)p
   - u(3,1,2)n >= -22
 MCLin9(3,1,2): - 22 y(1,3,2) + u(3,1,2)p - u(3,1,2)n <= -0
 MCLin10(3,1,2): 34 y(1,3,2) - v(3,1,2)p + v(3,1,2)n + u(3,1,2)p
   - u(3,1,2)n <= 34
 MCLinV(1,3,3): - x(1,3) + v(1,3,3)p - v(1,3,3)n + x(3,3) = -3
 MCLinV(3,1,3): x(1,3) + v(3,1,3)p - v(3,1,3)n - x(3,3) = -7
 MCLin1(1,3,3): u(1,3,3)p - u(1,3,3)n >= -0
 MCLin2(1,3,3): 33 y(3,1,3) + u(1,3,3)p - u(1,3,3)n >= -0
 MCLin3(1,3,3): - 23 y(3,1,3) - v(1,3,3)p + v(1,3,3)n + u(1,3,3)p
   - u(1,3,3)n >= -23
 MCLin4(1,3,3): - 23 y(3,1,3) + u(1,3,3)p - u(1,3,3)n <= -0
 MCLin5(1,3,3): 33 y(3,1,3) - v(1,3,3)p + v(1,3,3)n + u(1,3,3)p - u(1,3,3)n
   <= 33
 MCLin6(3,1,3): u(3,1,3)p - u(3,1,3)n >= -0
 MCLin7(3,1,3): 33 y(1,3,3) + u(3,1,3)p - u(3,1,3)n >= -0
 MCLin8(3,1,3): - 23 y(1,3,3) - v(3,1,3)p + v(3,1,3)n + u(3,1,3)p
   - u(3,1,3)n >= -23
 MCLin9(3,1,3): - 23 y(1,3,3) + u(3,1,3)p - u(3,1,3)n <= -0
 MCLin10(3,1,3): 33 y(1,3,3) - v(3,1,3)p + v(3,1,3)n + u(3,1,3)p
   - u(3,1,3)n <= 33
 MCLinV(2,3,0): - x(2,0) + v(2,3,0)p - v(2,3,0)n + x(3,0) = -9
 MCLinV(3,2,0): x(2,0) + v(3,2,0)p - v(3,2,0)n - x(3,0) = -1
 MCLin1(2,3,0): u(2,3,0)p - u(2,3,0)n >= -0
 MCLin2(2,3,0): 24 y(3,2,0) + u(2,3,0)p - u(2,3,0)n >= -0
 MCLin3(2,3,0): - 33 y(3,2,0) - v(2,3,0)p + v(2,3,0)n + u(2,3,0)p
   - u(2,3,0)n >= -33
 MCLin4(2,3,0): - 33 y(3,2,0) + u(2,3,0)p - u(2,3,0)n <= -0
 MCLin5(2,3,0): 24 y(3,2,0) - v(2,3,0)p + v(2,3,0)n + u(2,3,0)p - u(2,3,0)n
   <= 24
 MCLin6(3,2,0): u(3,2,0)p - u(3,2,0)n >= -0
 MCLin7(3,2,0): 43 y(2,3,0) + u(3,2,0)p - u(3,2,0)n >= -0
 MCLin8(3,2,0): - 14 y(2,3,0) - v(3,2,0)p + v(3,2,0)n + u(3,2,0)p
   - u(3,2,0)n >= -14
 MCLin9(3,2,0): - 14 y(2,3,0) + u(3,2,0)p - u(3,2,0)n <= -0
 MCLin10(3,2,0): 43 y(2,3,0) - v(3,2,0)p + v(3,2,0)n + u(3,2,0)p
   - u(3,2,0)n <= 43
 MCLinV(2,3,1): - x(2,1) + v(2,3,1)p - v(2,3,1)n + x(3,1) = -7
 MCLinV(3,2,1): x(2,1) + v(3,2,1)p - v(3,2,1)n - x(3,1) = -7
 MCLin1(2,3,1): u(2,3,1)p - u(2,3,1)n >= -0
 MCLin2(2,3,1): 19 y(3,2,1) + u(2,3,1)p - u(2,3,1)n >= -0
 MCLin3(2,3,1): - 38 y(3,2,1) - v(2,3,1)p + v(2,3,1)n + u(2,3,1)p
   - u(2,3,1)n >= -38
 MCLin4(2,3,1): - 38 y(3,2,1) + u(2,3,1)p - u(2,3,1)n <= -0
 MCLin5(2,3,1): 19 y(3,2,1) - v(2,3,1)p + v(2,3,1)n + u(2,3,1)p - u(2,3,1)n
   <= 19
 MCLin6(3,2,1): u(3,2,1)p - u(3,2,1)n >= -0
 MCLin7(3,2,1): 52 y(2,3,1) + u(3,2,1)p - u(3,2,1)n >= -0
 MCLin8(3,2,1): - 5 y(2,3,1) - v(3,2,1)p + v(3,2,1)n + u(3,2,1)p
   - u(3,2,1)n >= -5
 MCLin9(3,2,1): - 5 y(2,3,1) + u(3,2,1)p - u(3,2,1)n <= -0
 MCLin10(3,2,1): 52 y(2,3,1) - v(3,2,1)p + v(3,2,1)n + u(3,2,1)p
   - u(3,2,1)n <= 52
 MCLinV(2,3,2): - x(2,2) + v(2,3,2)p - v(2,3,2)n + x(3,2) = -10
 MCLinV(3,2,2): x(2,2) + v(3,2,2)p - v(3,2,2)n - x(3,2) = -3
 MCLin1(2,3,2): u(2,3,2)p - u(2,3,2)n >= -0
 MCLin2(2,3,2): 44 y(3,2,2) + u(2,3,2)p - u(2,3,2)n >= -0
 MCLin3(2,3,2): - 13 y(3,2,2) - v(2,3,2)p + v(2,3,2)n + u(2,3,2)p
   - u(2,3,2)n >= -13
 MCLin4(2,3,2): - 13 y(3,2,2) + u(2,3,2)p - u(2,3,2)n <= -0
 MCLin5(2,3,2): 44 y(3,2,2) - v(2,3,2)p + v(2,3,2)n + u(2,3,2)p - u(2,3,2)n
   <= 44
 MCLin6(3,2,2): u(3,2,2)p - u(3,2,2)n >= -0
 MCLin7(3,2,2): 26 y(2,3,2) + u(3,2,2)p - u(3,2,2)n >= -0
 MCLin8(3,2,2): - 31 y(2,3,2) - v(3,2,2)p + v(3,2,2)n + u(3,2,2)p
   - u(3,2,2)n >= -31
 MCLin9(3,2,2): - 31 y(2,3,2) + u(3,2,2)p - u(3,2,2)n <= -0
 MCLin10(3,2,2): 26 y(2,3,2) - v(3,2,2)p + v(3,2,2)n + u(3,2,2)p
   - u(3,2,2)n <= 26
 MCLinV(2,3,3): - x(2,3) + v(2,3,3)p - v(2,3,3)n + x(3,3) = -3
 MCLinV(3,2,3): x(2,3) + v(3,2,3)p - v(3,2,3)n - x(3,3) = -9
 MCLin1(2,3,3): u(2,3,3)p - u(2,3,3)n >= -0
 MCLin2(2,3,3): 34 y(3,2,3) + u(2,3,3)p - u(2,3,3)n >= -0
 MCLin3(2,3,3): - 23 y(3,2,3) - v(2,3,3)p + v(2,3,3)n + u(2,3,3)p
   - u(2,3,3)n >= -23
 MCLin4(2,3,3): - 23 y(3,2,3) + u(2,3,3)p - u(2,3,3)n <= -0
 MCLin5(2,3,3): 34 y(3,2,3) - v(2,3,3)p + v(2,3,3)n + u(2,3,3)p - u(2,3,3)n
   <= 34
 MCLin6(3,2,3): u(3,2,3)p - u(3,2,3)n >= -0
 MCLin7(3,2,3): 35 y(2,3,3) + u(3,2,3)p - u(3,2,3)n >= -0
 MCLin8(3,2,3): - 22 y(2,3,3) - v(3,2,3)p + v(3,2,3)n + u(3,2,3)p
   - u(3,2,3)n >= -22
 MCLin9(3,2,3): - 22 y(2,3,3) + u(3,2,3)p - u(3,2,3)n <= -0
 MCLin10(3,2,3): 35 y(2,3,3) - v(3,2,3)p + v(3,2,3)n + u(3,2,3)p
   - u(3,2,3)n <= 35
 makespan(0): C - x(0,1) >= 7
 makespan(1): C - x(1,2) >= 8
 makespan(2): C - x(2,0) >= 1
 makespan(3): C - x(3,2) >= 10
Bounds
 C <= 53
 x(0,0) <= 29
 v(0,1,0)p <= 28
 v(0,1,0)n <= 33
 v(1,0,0)p <= 22
 v(1,0,0)n <= 39
 u(0,1,0)p <= 530
 u(0,1,0)n <= 530
 u(1,0,0)p <= 530
 u(1,0,0)n <= 530
 v(0,2,0)p <= 9
 v(0,2,0)n <= 53
 v(2,0,0)p <= 42
 v(2,0,0)n <= 20
 u(0,2,0)p <= 530
 u(0,2,0)n <= 530
 u(2,0,0)p <= 530
 u(2,0,0)n <= 530
 v(0,3,0)p <= 10
 v(0,3,0)n <= 43
 v(3,0,0)p <= 24
 v(3,0,0)n <= 29
 u(0,3,0)p <= 530
 u(0,3,0)n <= 530
 u(3,0,0)p <= 530
 u(3,0,0)n <= 530
 x(1,0) <= 32
 v(1,2,0)p <= 12
 v(1,2,0)n <= 53
 v(2,1,0)p <= 51
 v(2,1,0)n <= 14
 u(1,2,0)p <= 530
 u(1,2,0)n <= 530
 u(2,1,0)p <= 530
 u(2,1,0)n <= 530
 v(1,3,0)p <= 13
 v(1,3,0)n <= 43
 v(3,1,0)p <= 33
 v(3,1,0)n <= 23
 u(1,3,0)p <= 530
 u(1,3,0)n <= 530
 u(3,1,0)p <= 530
 u(3,1,0)n <= 530
 19 <= x(2,0) <= 52
 v(2,3,0)p <= 33
 v(2,3,0)n <= 24
 v(3,2,0)p <= 14
 v(3,2,0)n <= 43
 u(2,3,0)p <= 530
 u(2,3,0)n <= 530
 u(3,2,0)p <= 530
 u(3,2,0)n <= 530
 10 <= x(3,0) <= 34
 17 <= x(0,1) <= 46
 v(0,1,1)p <= 33
 v(0,1,1)n <= 28
 v(1,0,1)p <= 16
 v(1,0,1)n <= 45
 u(0,1,1)p <= 530
 u(0,1,1)n <= 530
 u(1,0,1)p <= 530
 u(1,0,1)n <= 530
 v(0,2,1)p <= 27
 v(0,2,1)n <= 35
 v(2,0,1)p <= 21
 v(2,0,1)n <= 41
 u(0,2,1)p <= 530
 u(0,2,1)n <= 530
 u(2,0,1)p <= 530
 u(2,0,1)n <= 530
 v(0,3,1)p <= 39
 v(0,3,1)n <= 14
 v(3,0,1)p = 0
 v(3,0,1)n <= 53
 u(0,3,1)p <= 530
 u(0,3,1)n <= 530
 u(3,0,1)p <= 530
 u(3,0,1)n <= 530
 8 <= x(1,1) <= 40
 v(1,2,1)p <= 21
 v(1,2,1)n <= 44
 v(2,1,1)p <= 32
 v(2,1,1)n <= 33
 u(1,2,1)p <= 530
 u(1,2,1)n <= 530
 u(2,1,1)p <= 530
 u(2,1,1)n <= 530
 v(1,3,1)p <= 33
 v(1,3,1)n <= 23
 v(3,1,1)p <= 11
 v(3,1,1)n <= 45
 u(1,3,1)p <= 530
 u(1,3,1)n <= 530
 u(3,1,1)p <= 530
 u(3,1,1)n <= 530
 12 <= x(2,1) <= 45
 v(2,3,1)p <= 38
 v(2,3,1)n <= 19
 v(3,2,1)p <= 5
 v(3,2,1)n <= 52
 u(2,3,1)p <= 530
 u(2,3,1)n <= 530
 u(3,2,1)p <= 530
 u(3,2,1)n <= 530
 x(3,1) <= 24
 11 <= x(0,2) <= 40
 v(0,1,2)p <= 19
 v(0,1,2)n <= 42
 v(1,0,2)p <= 28
 v(1,0,2)n <= 33
 u(0,1,2)p <= 530
 u(0,1,2)n <= 530
 u(1,0,2)p <= 530
 u(1,0,2)n <= 530
 v(0,2,2)p <= 28
 v(0,2,2)n <= 34
 v(2,0,2)p <= 25
 v(2,0,2)n <= 37
 u(0,2,2)p <= 530
 u(0,2,2)n <= 530
 u(2,0,2)p <= 530
 u(2,0,2)n <= 530
 v(0,3,2)p <= 11
 v(0,3,2)n <= 42
 v(3,0,2)p <= 26
 v(3,0,2)n <= 27
 u(0,3,2)p <= 530
 u(0,3,2)n <= 530
 u(3,0,2)p <= 530
 u(3,0,2)n <= 530
 13 <= x(1,2) <= 45
 v(1,2,2)p <= 33
 v(1,2,2)n <= 32
 v(2,1,2)p <= 21
 v(2,1,2)n <= 44
 u(1,2,2)p <= 530
 u(1,2,2)n <= 530
 u(2,1,2)p <= 530
 u(2,1,2)n <= 530
 v(1,3,2)p <= 16
 v(1,3,2)n <= 40
 v(3,1,2)p <= 22
 v(3,1,2)n <= 34
 u(1,3,2)p <= 530
 u(1,3,2)n <= 530
 u(3,1,2)p <= 530
 u(3,1,2)n <= 530
 9 <= x(2,2) <= 42
 v(2,3,2)p <= 13
 v(2,3,2)n <= 44
 v(3,2,2)p <= 31
 v(3,2,2)n <= 26
 u(2,3,2)p <= 530
 u(2,3,2)n <= 530
 u(3,2,2)p <= 530
 u(3,2,2)n <= 530
 19 <= x(3,2) <= 43
 10 <= x(0,3) <= 39
 v(0,1,3)p <= 31
 v(0,1,3)n <= 30
 v(1,0,3)p <= 22
 v(1,0,3)n <= 39
 u(0,1,3)p <= 530
 u(0,1,3)n <= 530
 u(1,0,3)p <= 530
 u(1,0,3)n <= 530
 v(0,2,3)p <= 30
 v(0,2,3)n <= 32
 v(2,0,3)p <= 22
 v(2,0,3)n <= 40
 u(0,2,3)p <= 530
 u(0,2,3)n <= 530
 u(2,0,3)p <= 530
 u(2,0,3)n <= 530
 v(0,3,3)p <= 29
 v(0,3,3)n <= 24
 v(3,0,3)p <= 20
 v(3,0,3)n <= 33
 u(0,3,3)p <= 530
 u(0,3,3)n <= 530
 u(3,0,3)p <= 530
 u(3,0,3)n <= 530
 1 <= x(1,3) <= 33
 v(1,2,3)p <= 24
 v(1,2,3)n <= 41
 v(2,1,3)p <= 25
 v(2,1,3)n <= 40
 u(1,2,3)p <= 530
 u(1,2,3)n <= 530
 u(2,1,3)p <= 530
 u(2,1,3)n <= 530
 v(1,3,3)p <= 23
 v(1,3,3)n <= 33
 v(3,1,3)p <= 23
 v(3,1,3)n <= 33
 u(1,3,3)p <= 530
 u(1,3,3)n <= 530
 u(3,1,3)p <= 530
 u(3,1,3)n <= 530
 x(2,3) <= 33
 v(2,3,3)p <= 23
 v(2,3,3)n <= 34
 v(3,2,3)p <= 22
 v(3,2,3)n <= 35
 u(2,3,3)p <= 530
 u(2,3,3)n <= 530
 u(3,2,3)p <= 530
 u(3,2,3)n <= 530
 7 <= x(3,3) <= 31
Binaries
 y(0,1,0) y(1,0,0) y(0,2,0) y(2,0,0) y(0,3,0) y(3,0,0) y(1,2,0) y(2,1,0)
 y(1,3,0) y(3,1,0) y(2,3,0) y(3,2,0) y(0,1,1) y(1,0,1) y(0,2,1) y(2,0,1)
 y(0,3,1) y(3,0,1) y(1,2,1) y(2,1,1) y(1,3,1) y(3,1,1) y(2,3,1) y(3,2,1)
 y(0,1,2) y(1,0,2) y(0,2,2) y(2,0,2) y(0,3,2) y(3,0,2) y(1,2,2) y(2,1,2)
 y(1,3,2) y(3,1,2) y(2,3,2) y(3,2,2) y(0,1,3) y(1,0,3) y(0,2,3) y(2,0,3)
 y(0,3,3) y(3,0,3) y(1,2,3) y(2,1,3) y(1,3,3) y(3,1,3) y(2,3,3) y(3,2,3)
Generals
 C x(0,0) v(0,1,0)p v(0,1,0)n v(1,0,0)p v(1,0,0)n u(0,1,0)p u(0,1,0)n
 u(1,0,0)p u(1,0,0)n v(0,2,0)p v(0,2,0)n v(2,0,0)p v(2,0,0)n u(0,2,0)p
 u(0,2,0)n u(2,0,0)p u(2,0,0)n v(0,3,0)p v(0,3,0)n v(3,0,0)p v(3,0,0)n
 u(0,3,0)p u(0,3,0)n u(3,0,0)p u(3,0,0)n x(1,0) v(1,2,0)p v(1,2,0)n
 v(2,1,0)p v(2,1,0)n u(1,2,0)p u(1,2,0)n u(2,1,0)p u(2,1,0)n v(1,3,0)p
 v(1,3,0)n v(3,1,0)p v(3,1,0)n u(1,3,0)p u(1,3,0)n u(3,1,0)p u(3,1,0)n
 x(2,0) v(2,3,0)p v(2,3,0)n v(3,2,0)p v(3,2,0)n u(2,3,0)p u(2,3,0)n
 u(3,2,0)p u(3,2,0)n x(3,0) x(0,1) v(0,1,1)p v(0,1,1)n v(1,0,1)p v(1,0,1)n
 u(0,1,1)p u(0,1,1)n u(1,0,1)p u(1,0,1)n v(0,2,1)p v(0,2,1)n v(2,0,1)p
 v(2,0,1)n u(0,2,1)p u(0,2,1)n u(2,0,1)p u(2,0,1)n v(0,3,1)p v(0,3,1)n
 v(3,0,1)p v(3,0,1)n u(0,3,1)p u(0,3,1)n u(3,0,1)p u(3,0,1)n x(1,1)
 v(1,2,1)p v(1,2,1)n v(2,1,1)p v(2,1,1)n u(1,2,1)p u(1,2,1)n u(2,1,1)p
 u(2,1,1)n v(1,3,1)p v(1,3,1)n v(3,1,1)p v(3,1,1)n u(1,3,1)p u(1,3,1)n
 u(3,1,1)p u(3,1,1)n x(2,1) v(2,3,1)p v(2,3,1)n v(3,2,1)p v(3,2,1)n
 u(2,3,1)p u(2,3,1)n u(3,2,1)p u(3,2,1)n x(3,1) x(0,2) v(0,1,2)p v(0,1,2)n
 v(1,0,2)p v(1,0,2)n u(0,1,2)p u(0,1,2)n u(1,0,2)p u(1,0,2)n v(0,2,2)p
 v(0,2,2)n v(2,0,2)p v(2,0,2)n u(0,2,2)p u(0,2,2)n u(2,0,2)p u(2,0,2)n
 v(0,3,2)p v(0,3,2)n v(3,0,2)p v(3,0,2)n u(0,3,2)p u(0,3,2)n u(3,0,2)p
 u(3,0,2)n x(1,2) v(1,2,2)p v(1,2,2)n v(2,1,2)p v(2,1,2)n u(1,2,2)p
 u(1,2,2)n u(2,1,2)p u(2,1,2)n v(1,3,2)p v(1,3,2)n v(3,1,2)p v(3,1,2)n
 u(1,3,2)p u(1,3,2)n u(3,1,2)p u(3,1,2)n x(2,2) v(2,3,2)p v(2,3,2)n
 v(3,2,2)p v(3,2,2)n u(2,3,2)p u(2,3,2)n u(3,2,2)p u(3,2,2)n x(3,2) x(0,3)
 v(0,1,3)p v(0,1,3)n v(1,0,3)p v(1,0,3)n u(0,1,3)p u(0,1,3)n u(1,0,3)p
 u(1,0,3)n v(0,2,3)p v(0,2,3)n v(2,0,3)p v(2,0,3)n u(0,2,3)p u(0,2,3)n
 u(2,0,3)p u(2,0,3)n v(0,3,3)p v(0,3,3)n v(3,0,3)p v(3,0,3)n u(0,3,3)p
 u(0,3,3)n u(3,0,3)p u(3,0,3)n x(1,3) v(1,2,3)p v(1,2,3)n v(2,1,3)p
 v(2,1,3)n u(1,2,3)p u(1,2,3)n u(2,1,3)p u(2,1,3)n v(1,3,3)p v(1,3,3)n
 v(3,1,3)p v(3,1,3)n u(1,3,3)p u(1,3,3)n u(3,1,3)p u(3,1,3)n x(2,3)
 v(2,3,3)p v(2,3,3)n v(3,2,3)p v(3,2,3)n u(2,3,3)p u(2,3,3)n u(3,2,3)p
 u(3,2,3)n x(3,3)
End
