\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n - 2 x1p + 2 x1n + x3p - x3n + 4 x4p - 4 x4n
Subject To
 constr(0): - 5 x0p + 5 x0n - 2 x1p + 2 x1n - x4p + x4n = 3
 constr(1): 5 x0p - 5 x0n - 2 x1p + 2 x1n + 5 x2p - 5 x2n - 3 x3p + 3 x3n
   = 13
 constr(2): - x1p + x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n - 4 x4p + 4 x4n
   >= -8
Bounds
 x0p = 0
 x0n <= 3
 x1p <= 4
 x1n <= 5
 x2p <= 2
 x2n <= 2
 x3p = 0
 x3n <= 5
 x4p <= 1
 x4n <= 2
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4p x4n
End
