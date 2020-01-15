\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  0.8 a(x1) + 2 a(x3) + 5 a(x2) - a0
Subject To
 f_x1: - u(constr(1)) - a(x1) - f(x1) = -0
 f_x3: - 2 u(constr(0)) - a(x3) - f(x3) = -0
 f_x2: 3 u(constr(0)) + 3 u(constr(1)) - a(x2) - f(x2) = -0
 f: 11 u(constr(0)) + 14.2 u(constr(1)) - a0 - f0 = -0
Bounds
 u(constr(0)) <= 0.99
 u(constr(1)) <= 0.99
 a(x1) >= -100000
 a(x3) >= -100000
 a(x2) >= -100000
 f(x1) <= 0.99
 f(x3) <= 0.99
 f(x2) <= 0.99
 a0 >= -100000
 f0 <= 0.99
Generals
 a(x1) a(x3) a(x2) a0
End
