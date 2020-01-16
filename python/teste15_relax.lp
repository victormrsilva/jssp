\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + x1p - x1n - 2 x3p + 2 x3n
Subject To
 constr(0): x1p - x1n + 3 x2p - 3 x2n + 2 x3p - 2 x3n <= 3
 constr(1): - 3 x0p + 3 x0n - 2 x1p + 2 x1n - x2p + x2n = 6
Bounds
 x0p <= 5
 x0n <= 3
 x1p <= 5
 x1n <= 2
 x2p <= 2
 x2n <= 5
 x3p = 0
 x3n <= 2
End
