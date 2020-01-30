\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x1p + 2 x1n + 4 x4p - 4 x4n
Subject To
 constr(0): 3 x0p - 3 x0n + x1p - x1n + 3 x3p - 3 x3n - x4p + x4n >= -19
 constr(1): 5 x0p - 5 x0n - 4 x2p + 4 x2n + 3 x4p - 3 x4n >= -0
 constr(2): 2 x0p - 2 x0n + 4 x2p - 4 x2n - 5 x3p + 5 x3n + 3 x4p - 3 x4n
   = 23
Bounds
 x0p <= 4
 x0n <= 3
 x1p <= 1
 x1n <= 3
 x2p <= 3
 x2n <= 2
 x3p <= 4
 x3n <= 5
 x4p <= 5
 x4n <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4p x4n
End
