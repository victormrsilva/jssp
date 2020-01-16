\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x1p + 3 x1n - 5 x2p + 5 x2n + 5 x3p - 5 x3n
Subject To
 constr(0): - x0p + x0n + 5 x1p - 5 x1n + x2p - x2n <= 6
 constr(1): - x0p + x0n - x1p + x1n - 2 x2p + 2 x2n + 4 x3p - 4 x3n >= -20
Bounds
 x0p <= 1
 x0n <= 2
 x1p <= 3
 x1n <= 4
 x2p = 0
 x2n <= 2
 x3p <= 3
 x3n <= 3
End
