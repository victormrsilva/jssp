\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  0.25 a(x4n) - a0
Subject To
 f_x4n: 5 u(constr(0)) - 4 u(constr(1)) - a(x4n) - f(x4n) = -0
 f: 64 u(constr(0)) - u(constr(1)) - a0 - f0 = -0
Bounds
 u(constr(0)) <= 0.99
 u(constr(1)) <= 0.99
 a(x4n) >= -100000
 f(x4n) <= 0.99
 a0 >= -100000
 f0 <= 0.99
Generals
 a(x4n) a0
End
