\Problem name: compact

Minimize
OBJROW: C
Subject To
ord(0,1):  - x(0,1) + x(0,2) >= 7
ord(0,2):  x(0,0) - x(0,2) >= 7
ord(1,1):  x(1,0) - x(1,1) >= 4
ord(1,2):  - x(1,0) + x(1,2) >= 3
ord(2,1):  - x(2,1) + x(2,2) >= 7
ord(2,2):  x(2,0) - x(2,2) >= 6
triangle2(0,1,0):  y(0,1,0) + y(1,0,0) = 1
triangle2(0,1,1):  y(0,1,1) + y(1,0,1) = 1
triangle2(0,1,2):  y(0,1,2) + y(1,0,2) = 1
triangle2(0,2,0):  y(0,2,0) + y(2,0,0) = 1
triangle2(0,2,1):  y(0,2,1) + y(2,0,1) = 1
triangle2(0,2,2):  y(0,2,2) + y(2,0,2) = 1
triangle2(1,2,0):  y(1,2,0) + y(2,1,0) = 1
triangle2(1,2,1):  y(1,2,1) + y(2,1,1) = 1
triangle2(1,2,2):  y(1,2,2) + y(2,1,2) = 1
MCLinV(0,1,0):  x(0,0) - v(0,1,0) - x(1,0) = 3
MCLinV(1,0,0):  - x(0,0) - v(1,0,0) + x(1,0) = 1
MCLin1(0,1,0):  u(0,1,0) >= 0
MCLin2(0,1,0):  -23 y(1,0,0) - u(0,1,0) <= 0
MCLin3(0,1,0):  36 y(1,0,0) + v(0,1,0) - u(0,1,0) <= 36
MCLin4(0,1,0):  36 y(1,0,0) - u(0,1,0) >= 0
MCLin5(0,1,0):  -23 y(1,0,0) + v(0,1,0) - u(0,1,0) >= -23
MCLin6(1,0,0):  u(1,0,0) >= 0
MCLin7(1,0,0):  -40 y(0,1,0) - u(1,0,0) <= 0
MCLin8(1,0,0):  19 y(0,1,0) + v(1,0,0) - u(1,0,0) <= 19
MCLin9(1,0,0):  19 y(0,1,0) - u(1,0,0) >= 0
MCLin10(1,0,0):  -40 y(0,1,0) + v(1,0,0) - u(1,0,0) >= -40
MCLinV(0,1,1):  x(0,1) - v(0,1,1) - x(1,1) = 4
MCLinV(1,0,1):  - x(0,1) - v(1,0,1) + x(1,1) = 7
MCLin1(0,1,1):  u(0,1,1) >= 0
MCLin2(0,1,1):  -34 y(1,0,1) - u(0,1,1) <= 0
MCLin3(0,1,1):  25 y(1,0,1) + v(0,1,1) - u(0,1,1) <= 25
MCLin4(0,1,1):  25 y(1,0,1) - u(0,1,1) >= 0
MCLin5(0,1,1):  -34 y(1,0,1) + v(0,1,1) - u(0,1,1) >= -34
MCLin6(1,0,1):  u(1,0,1) >= 0
MCLin7(1,0,1):  -36 y(0,1,1) - u(1,0,1) <= 0
MCLin8(1,0,1):  23 y(0,1,1) + v(1,0,1) - u(1,0,1) <= 23
MCLin9(1,0,1):  23 y(0,1,1) - u(1,0,1) >= 0
MCLin10(1,0,1):  -36 y(0,1,1) + v(1,0,1) - u(1,0,1) >= -36
MCLinV(0,1,2):  x(0,2) - v(0,1,2) - x(1,2) = 7
MCLinV(1,0,2):  - x(0,2) - v(1,0,2) + x(1,2) = 7
MCLin1(0,1,2):  u(0,1,2) >= 0
MCLin2(0,1,2):  -37 y(1,0,2) - u(0,1,2) <= 0
MCLin3(0,1,2):  22 y(1,0,2) + v(0,1,2) - u(0,1,2) <= 22
MCLin4(0,1,2):  22 y(1,0,2) - u(0,1,2) >= 0
MCLin5(0,1,2):  -37 y(1,0,2) + v(0,1,2) - u(0,1,2) >= -37
MCLin6(1,0,2):  u(1,0,2) >= 0
MCLin7(1,0,2):  -36 y(0,1,2) - u(1,0,2) <= 0
MCLin8(1,0,2):  23 y(0,1,2) + v(1,0,2) - u(1,0,2) <= 23
MCLin9(1,0,2):  23 y(0,1,2) - u(1,0,2) >= 0
MCLin10(1,0,2):  -36 y(0,1,2) + v(1,0,2) - u(1,0,2) >= -36
MCLinV(0,2,0):  x(0,0) - v(0,2,0) - x(2,0) = 2
MCLinV(2,0,0):  - x(0,0) - v(2,0,0) + x(2,0) = 1
MCLin1(0,2,0):  u(0,2,0) >= 0
MCLin2(0,2,0):  -30 y(2,0,0) - u(0,2,0) <= 0
MCLin3(0,2,0):  28 y(2,0,0) + v(0,2,0) - u(0,2,0) <= 28
MCLin4(0,2,0):  28 y(2,0,0) - u(0,2,0) >= 0
MCLin5(0,2,0):  -30 y(2,0,0) + v(0,2,0) - u(0,2,0) >= -30
MCLin6(2,0,0):  u(2,0,0) >= 0
MCLin7(2,0,0):  -31 y(0,2,0) - u(2,0,0) <= 0
MCLin8(2,0,0):  27 y(0,2,0) + v(2,0,0) - u(2,0,0) <= 27
MCLin9(2,0,0):  27 y(0,2,0) - u(2,0,0) >= 0
MCLin10(2,0,0):  -31 y(0,2,0) + v(2,0,0) - u(2,0,0) >= -31
MCLinV(0,2,1):  x(0,1) - v(0,2,1) - x(2,1) = 7
MCLinV(2,0,1):  - x(0,1) - v(2,0,1) + x(2,1) = 7
MCLin1(0,2,1):  u(0,2,1) >= 0
MCLin2(0,2,1):  -36 y(2,0,1) - u(0,2,1) <= 0
MCLin3(0,2,1):  22 y(2,0,1) + v(0,2,1) - u(0,2,1) <= 22
MCLin4(0,2,1):  22 y(2,0,1) - u(0,2,1) >= 0
MCLin5(0,2,1):  -36 y(2,0,1) + v(0,2,1) - u(0,2,1) >= -36
MCLin6(2,0,1):  u(2,0,1) >= 0
MCLin7(2,0,1):  -36 y(0,2,1) - u(2,0,1) <= 0
MCLin8(2,0,1):  22 y(0,2,1) + v(2,0,1) - u(2,0,1) <= 22
MCLin9(2,0,1):  22 y(0,2,1) - u(2,0,1) >= 0
MCLin10(2,0,1):  -36 y(0,2,1) + v(2,0,1) - u(2,0,1) >= -36
MCLinV(0,2,2):  x(0,2) - v(0,2,2) - x(2,2) = 6
MCLinV(2,0,2):  - x(0,2) - v(2,0,2) + x(2,2) = 7
MCLin1(0,2,2):  u(0,2,2) >= 0
MCLin2(0,2,2):  -35 y(2,0,2) - u(0,2,2) <= 0
MCLin3(0,2,2):  23 y(2,0,2) + v(0,2,2) - u(0,2,2) <= 23
MCLin4(0,2,2):  23 y(2,0,2) - u(0,2,2) >= 0
MCLin5(0,2,2):  -35 y(2,0,2) + v(0,2,2) - u(0,2,2) >= -35
MCLin6(2,0,2):  u(2,0,2) >= 0
MCLin7(2,0,2):  -36 y(0,2,2) - u(2,0,2) <= 0
MCLin8(2,0,2):  22 y(0,2,2) + v(2,0,2) - u(2,0,2) <= 22
MCLin9(2,0,2):  22 y(0,2,2) - u(2,0,2) >= 0
MCLin10(2,0,2):  -36 y(0,2,2) + v(2,0,2) - u(2,0,2) >= -36
MCLinV(1,2,0):  x(1,0) - v(1,2,0) - x(2,0) = 2
MCLinV(2,1,0):  - x(1,0) - v(2,1,0) + x(2,0) = 3
MCLin1(1,2,0):  u(1,2,0) >= 0
MCLin2(1,2,0):  -40 y(2,1,0) - u(1,2,0) <= 0
MCLin3(1,2,0):  19 y(2,1,0) + v(1,2,0) - u(1,2,0) <= 19
MCLin4(1,2,0):  19 y(2,1,0) - u(1,2,0) >= 0
MCLin5(1,2,0):  -40 y(2,1,0) + v(1,2,0) - u(1,2,0) >= -40
MCLin6(2,1,0):  u(2,1,0) >= 0
MCLin7(2,1,0):  -24 y(1,2,0) - u(2,1,0) <= 0
MCLin8(2,1,0):  35 y(1,2,0) + v(2,1,0) - u(2,1,0) <= 35
MCLin9(2,1,0):  35 y(1,2,0) - u(2,1,0) >= 0
MCLin10(2,1,0):  -24 y(1,2,0) + v(2,1,0) - u(2,1,0) >= -24
MCLinV(1,2,1):  x(1,1) - v(1,2,1) - x(2,1) = 7
MCLinV(2,1,1):  - x(1,1) - v(2,1,1) + x(2,1) = 4
MCLin1(1,2,1):  u(1,2,1) >= 0
MCLin2(1,2,1):  -36 y(2,1,1) - u(1,2,1) <= 0
MCLin3(1,2,1):  23 y(2,1,1) + v(1,2,1) - u(1,2,1) <= 23
MCLin4(1,2,1):  23 y(2,1,1) - u(1,2,1) >= 0
MCLin5(1,2,1):  -36 y(2,1,1) + v(1,2,1) - u(1,2,1) >= -36
MCLin6(2,1,1):  u(2,1,1) >= 0
MCLin7(2,1,1):  -34 y(1,2,1) - u(2,1,1) <= 0
MCLin8(2,1,1):  25 y(1,2,1) + v(2,1,1) - u(2,1,1) <= 25
MCLin9(2,1,1):  25 y(1,2,1) - u(2,1,1) >= 0
MCLin10(2,1,1):  -34 y(1,2,1) + v(2,1,1) - u(2,1,1) >= -34
MCLinV(1,2,2):  x(1,2) - v(1,2,2) - x(2,2) = 6
MCLinV(2,1,2):  - x(1,2) - v(2,1,2) + x(2,2) = 7
MCLin1(1,2,2):  u(1,2,2) >= 0
MCLin2(1,2,2):  -35 y(2,1,2) - u(1,2,2) <= 0
MCLin3(1,2,2):  24 y(2,1,2) + v(1,2,2) - u(1,2,2) <= 24
MCLin4(1,2,2):  24 y(2,1,2) - u(1,2,2) >= 0
MCLin5(1,2,2):  -35 y(2,1,2) + v(1,2,2) - u(1,2,2) >= -35
MCLin6(2,1,2):  u(2,1,2) >= 0
MCLin7(2,1,2):  -37 y(1,2,2) - u(2,1,2) <= 0
MCLin8(2,1,2):  22 y(1,2,2) + v(2,1,2) - u(2,1,2) <= 22
MCLin9(2,1,2):  22 y(1,2,2) - u(2,1,2) >= 0
MCLin10(2,1,2):  -37 y(1,2,2) + v(2,1,2) - u(2,1,2) >= -37
makespan(0):  C - x(0,0) >= 1
makespan(1):  C - x(1,2) >= 7
makespan(2):  C - x(2,0) >= 2
cut_clique_1(012,0):  0.04651 x(0,0) + 0.02326 x(2,0) >= 1
cut_clique_1(012,1):  0.25000 x(0,1) + 0.14286 x(1,1) >= 1
cut_clique_1(012,2):  0.02564 x(0,2) + 0.02564 x(1,2) + 0.02198 x(2,2) >= 1
cut_basic_best(02,0):  x(0,0) + 2 x(2,0) >= 41
cut_two_jobs(0,2,0):  2 x(0,0) + x(2,0) >= 43
cut_basic_best(012,1):  7 x(0,1) + 4 x(1,1) + 7 x(2,1) >= 105
cut_two_jobs(0,2,1):  7 x(0,1) + 7 x(2,1) >= 49
cut_basic_plus_epsilon_best(12,1,0):  4 x(1,1) + 7 x(2,1) >= 28
cut_basic_plus_epsilon_best(02,1,1):  7 x(0,1) + 7 x(2,1) >= 49
cut_basic_plus_epsilon_best(01,1,2):  7 x(0,1) + 4 x(1,1) >= 28
cut_basic_best(012,2):  7 x(0,2) + 7 x(1,2) + 6 x(2,2) >= 273
cut_two_jobs(0,1,2):  7 x(0,2) + 7 x(1,2) >= 147
cut_basic_plus_epsilon_best(12,2,0):  7 x(1,2) + 6 x(2,2) >= 133
cut_basic_plus_epsilon_best(02,2,1):  7 x(0,2) + 6 x(2,2) >= 133
cut_basic_plus_epsilon_best(01,2,2):  7 x(0,2) + 7 x(1,2) >= 147
Bounds
 0 <= C <= 44
 14 <= x(0,0) <= 43
 0 <= y(0,1,0) <= 1
 0 <= y(1,0,0) <= 1
 -23 <= v(0,1,0) <= 36
 -40 <= v(1,0,0) <= 19
 -440 <= u(0,1,0)
 -440 <= u(1,0,0)
 0 <= y(0,2,0) <= 1
 0 <= y(2,0,0) <= 1
 -30 <= v(0,2,0) <= 28
 -31 <= v(2,0,0) <= 27
 -440 <= u(0,2,0)
 -440 <= u(2,0,0)
 4 <= x(1,0) <= 34
 0 <= y(1,2,0) <= 1
 0 <= y(2,1,0) <= 1
 -40 <= v(1,2,0) <= 19
 -24 <= v(2,1,0) <= 35
 -440 <= u(1,2,0)
 -440 <= u(2,1,0)
 13 <= x(2,0) <= 42
 0 <= x(0,1) <= 29
 0 <= y(0,1,1) <= 1
 0 <= y(1,0,1) <= 1
 -34 <= v(0,1,1) <= 25
 -36 <= v(1,0,1) <= 23
 -440 <= u(0,1,1)
 -440 <= u(1,0,1)
 0 <= y(0,2,1) <= 1
 0 <= y(2,0,1) <= 1
 -36 <= v(0,2,1) <= 22
 -36 <= v(2,0,1) <= 22
 -440 <= u(0,2,1)
 -440 <= u(2,0,1)
 0 <= x(1,1) <= 30
 0 <= y(1,2,1) <= 1
 0 <= y(2,1,1) <= 1
 -36 <= v(1,2,1) <= 23
 -34 <= v(2,1,1) <= 25
 -440 <= u(1,2,1)
 -440 <= u(2,1,1)
 0 <= x(2,1) <= 29
 7 <= x(0,2) <= 36
 0 <= y(0,1,2) <= 1
 0 <= y(1,0,2) <= 1
 -37 <= v(0,1,2) <= 22
 -36 <= v(1,0,2) <= 23
 -440 <= u(0,1,2)
 -440 <= u(1,0,2)
 0 <= y(0,2,2) <= 1
 0 <= y(2,0,2) <= 1
 -35 <= v(0,2,2) <= 23
 -36 <= v(2,0,2) <= 22
 -440 <= u(0,2,2)
 -440 <= u(2,0,2)
 7 <= x(1,2) <= 37
 0 <= y(1,2,2) <= 1
 0 <= y(2,1,2) <= 1
 -35 <= v(1,2,2) <= 24
 -37 <= v(2,1,2) <= 22
 -440 <= u(1,2,2)
 -440 <= u(2,1,2)
 7 <= x(2,2) <= 36
End
