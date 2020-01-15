\ Model chvatal
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  1.333333333333333 a(x0n) + 2 a(x3n) + 2.666666666666666 a(x4p) + 3 a(x1p)
   - a0
Subject To
 f_x2: - 4 u(constr(0)) - 2 u(constr(1)) - 5 u(constr(2)) - a(x2) - f(x2)
   = -0
 f_x0p: 3 u(constr(0)) + 3 u(constr(2)) + u(constr(3)) - a(x0p) - f(x0p)
   = -0
 f_x0n: - 3 u(constr(0)) - 3 u(constr(2)) - u(constr(3)) - a(x0n) - f(x0n)
   = -0
 f_x3p: u(constr(0)) + 5 u(constr(1)) + 2 u(constr(3)) - a(x3p) - f(x3p)
   = -0
 f_x3n: - u(constr(0)) - 5 u(constr(1)) - 2 u(constr(3)) - a(x3n) - f(x3n)
   = -0
 f_x4p: 3 u(constr(0)) - 5 u(constr(1)) - 4 u(constr(2)) + 2 u(constr(3))
   - a(x4p) - f(x4p) = -0
 f_x4n: - 3 u(constr(0)) + 5 u(constr(1)) + 4 u(constr(2)) - 2 u(constr(3))
   - a(x4n) - f(x4n) = -0
 f_x1p: - 5 u(constr(2)) - a(x1p) - f(x1p) = -0
 f_x1n: 5 u(constr(2)) - a(x1n) - f(x1n) = -0
 f: 2 u(constr(0)) + u(constr(1)) + 18 u(constr(2)) - a0 - f0 = -0
Bounds
 u(constr(0)) <= 0.99
 u(constr(1)) <= 0.99
 u(constr(2)) <= 0.99
 u(constr(3)) <= 0.99
 a(x0p) >= -100000
 a(x0n) >= -100000
 f(x0p) <= 0.99
 f(x0n) <= 0.99
 a(x2) >= -100000
 f(x2) <= 0.99
 a(x3p) >= -100000
 a(x3n) >= -100000
 f(x3p) <= 0.99
 f(x3n) <= 0.99
 a(x4p) >= -100000
 a(x4n) >= -100000
 f(x4p) <= 0.99
 f(x4n) <= 0.99
 a(x1p) >= -100000
 a(x1n) >= -100000
 f(x1p) <= 0.99
 f(x1n) <= 0.99
 a0 >= -100000
 f0 <= 0.99
Generals
 a(x0p) a(x0n) a(x2) a(x3p) a(x3n) a(x4p) a(x4n) a(x1p) a(x1n) a0
End
