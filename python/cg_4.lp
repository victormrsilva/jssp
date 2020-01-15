\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  a(x0n) + 1.846153846153846 a(x1p) + 0.6153846153846156 a(x1n) + 5 a(x2p)
   + 1.538461538461538 a(x4p) + 5 a(x4n) - a0
Subject To
 f_x0p: - u(constr(0)) + 4 u(constr(1)) + 5 u(constr(2)) - a(x0p) - f(x0p)
   = -0
 f_x0n: u(constr(0)) - 4 u(constr(1)) - 5 u(constr(2)) - a(x0n) - f(x0n)
   = -0
 f_x1p: 3 u(constr(0)) - 5 u(constr(1)) - 4 u(constr(2)) - a(x1p) - f(x1p)
   = -0
 f_x1n: - 3 u(constr(0)) + 5 u(constr(1)) + 4 u(constr(2)) - a(x1n)
   - f(x1n) = -0
 f_x2p: - 4 u(constr(0)) - 5 u(constr(1)) + 2 u(constr(2)) - a(x2p)
   - f(x2p) = -0
 f_x2n: 4 u(constr(0)) + 5 u(constr(1)) - 2 u(constr(2)) - a(x2n) - f(x2n)
   = -0
 f_x3p: - 2 u(constr(0)) - 2 u(constr(2)) - a(x3p) - f(x3p) = -0
 f_x3n: 2 u(constr(0)) + 2 u(constr(2)) - a(x3n) - f(x3n) = -0
 f_x4p: - 5 u(constr(0)) + 4 u(constr(1)) - 2 u(constr(2)) - a(x4p)
   - f(x4p) = -0
 f_x4n: 5 u(constr(0)) - 4 u(constr(1)) + 2 u(constr(2)) - a(x4n) - f(x4n)
   = -0
 f: 2 u(constr(0)) + 16 u(constr(1)) + 7 u(constr(2)) - a0 - f0 = -0
Bounds
 u(constr(0)) <= 0.99
 u(constr(1)) <= 0.99
 u(constr(2)) <= 0.99
 a(x0p) >= -100000
 a(x0n) >= -100000
 a(x1p) >= -100000
 a(x1n) >= -100000
 a(x2p) >= -100000
 a(x2n) >= -100000
 a(x3p) >= -100000
 a(x3n) >= -100000
 a(x4p) >= -100000
 a(x4n) >= -100000
 f(x0p) <= 0.99
 f(x0n) <= 0.99
 f(x1p) <= 0.99
 f(x1n) <= 0.99
 f(x2p) <= 0.99
 f(x2n) <= 0.99
 f(x3p) <= 0.99
 f(x3n) <= 0.99
 f(x4p) <= 0.99
 f(x4n) <= 0.99
 a0 >= -100000
 f0 <= 0.99
Generals
 a(x0p) a(x0n) a(x1p) a(x1n) a(x2p) a(x2n) a(x3p) a(x3n) a(x4p) a(x4n) a0
End
