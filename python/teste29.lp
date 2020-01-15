\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n - x1p + x1n + x2p - x2n - 2 x3p + 2 x3n
Subject To
 constr(0): - 2 x1p + 2 x1n - 2 x2p + 2 x2n + 5 x3p - 5 x3n = 9
 constr(1): 4 x0p - 4 x0n - 5 x1p + 5 x1n + 2 x2p - 2 x2n - x3p + x3n
   >= -0
Bounds
 x0p <= 4
 x0n <= 1
 x1p <= 4
 x1n <= 4
 x2p = 0
 x2n <= 1
 x3p <= 3
 x3n <= 1
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
