\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n - 2 x1p + 2 x1n + 5 x2p - 5 x2n + 5 x3p - 5 x3n
Subject To
 constr(0): 3 x1p - 3 x1n - x2p + x2n - 4 x3p + 4 x3n <= 15
 constr(1): - 5 x0p + 5 x0n + 4 x1p - 4 x1n + 4 x2p - 4 x2n + 3 x3p - 3 x3n
   >= -0
Bounds
 x0p = 0
 x0n <= 3
 x1p <= 4
 x1n <= 2
 x2p <= 5
 x2n <= 2
 x3p <= 4
 x3n <= 2
End
