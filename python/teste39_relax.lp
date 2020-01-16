\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n + 3 x1p - 3 x1n - 4 x2p + 4 x2n
Subject To
 constr(0): 4 x0p - 4 x0n + 3 x1p - 3 x1n - 3 x2p + 3 x2n = 2
 constr(1): - x0p + x0n - 4 x1p + 4 x1n >= -0
Bounds
 x0p <= 5
 x0n <= 5
 x1p <= 1
 x1n <= 1
 x2p = 0
 x2n <= 4
End
