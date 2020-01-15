\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  4 a(x0n) + 3 a(x1n) + 5 a(x2p) + 0.25 a(x4n) - a0
Subject To
 f_x0p: 5 u(constr(0)) - a(x0p) - f(x0p) = -0
 f_x0n: - 5 u(constr(0)) - a(x0n) - f(x0n) = -0
 f_x1p: 5 u(constr(0)) + 2 u(constr(1)) - a(x1p) - f(x1p) = -0
 f_x1n: - 5 u(constr(0)) - 2 u(constr(1)) - a(x1n) - f(x1n) = -0
 f_x2p: - 3 u(constr(0)) + 3 u(constr(1)) - a(x2p) - f(x2p) = -0
 f_x2n: 3 u(constr(0)) - 3 u(constr(1)) - a(x2n) - f(x2n) = -0
 f_x3: - 5 u(constr(0)) + u(constr(1)) - a(x3) - f(x3) = -0
 f_x4p: - 5 u(constr(0)) + 4 u(constr(1)) - a(x4p) - f(x4p) = -0
 f_x4n: 5 u(constr(0)) - 4 u(constr(1)) - a(x4n) - f(x4n) = -0
 f: 14 u(constr(0)) + 8 u(constr(1)) - a0 - f0 = -0
Bounds
 u(constr(0)) <= 0.99
 u(constr(1)) <= 0.99
 a(x0p) >= -100000
 f(x0p) <= 0.99
 a(x0n) >= -100000
 f(x0n) <= 0.99
 a(x1p) >= -100000
 f(x1p) <= 0.99
 a(x1n) >= -100000
 f(x1n) <= 0.99
 a(x2p) >= -100000
 f(x2p) <= 0.99
 a(x2n) >= -100000
 f(x2n) <= 0.99
 a(x3) >= -100000
 f(x3) <= 0.99
 a(x4p) >= -100000
 f(x4p) <= 0.99
 a(x4n) >= -100000
 f(x4n) <= 0.99
 a0 >= -100000
 f0 <= 0.99
Generals
 a(x0p) a(x0n) a(x1p) a(x1n) a(x2p) a(x2n) a(x3) a(x4p) a(x4n) a0
End
