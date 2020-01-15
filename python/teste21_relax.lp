\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n + 3 x1p - 3 x1n - x2p + x2n
Subject To
 constr(0): 2 x0p - 2 x0n + x1p - x1n + 3 x2p - 3 x2n = -0
 constr(1): - 4 x1p + 4 x1n - 5 x2p + 5 x2n <= 4
Bounds
 x0p <= 5
 x0n <= 4
 x1p <= 2
 x1n <= 4
 x2p <= 2
 x2n <= 2
End
