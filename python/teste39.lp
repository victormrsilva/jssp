\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n - x3p + x3n
Subject To
 constr(0): 3 x0p - 3 x0n + x2p - x2n = 3
 constr(1): 3 x0p - 3 x0n + 5 x2p - 5 x2n - 3 x3p + 3 x3n >= -0
 constr(2): - 3 x0p + 3 x0n + x3p - x3n <= 13
 constr(3): 4 x0p - 4 x0n + 4 x1p - 4 x1n + 5 x2p - 5 x2n - 5 x3p + 5 x3n
   = 18
 constr(4): x1p - x1n - 4 x2p + 4 x2n - 4 x3p + 4 x3n <= 18
Bounds
 x0p <= 1
 x0n <= 3
 x1p <= 5
 x1n <= 1
 x2p <= 2
 x2n <= 4
 x3p <= 5
 x3n <= 2
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
