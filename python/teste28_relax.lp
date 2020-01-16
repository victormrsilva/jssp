\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x1p + 3 x1n + 3 x2p - 3 x2n + x3p - x3n
Subject To
 constr(0): 2 x0p - 2 x0n + 2 x1p - 2 x1n - 2 x3p + 2 x3n <= 7
 constr(1): 3 x1p - 3 x1n - x3p + x3n >= -0
 constr(2): x0p - x0n - 2 x2p + 2 x2n + 2 x3p - 2 x3n = 3
 constr(3): 3 x0p - 3 x0n + x1p - x1n - x2p + x2n - 3 x3p + 3 x3n <= 16
Bounds
 x0p <= 5
 x0n <= 4
 x1p <= 3
 x1n <= 1
 x2p <= 5
 x2n <= 2
 x3p <= 5
 x3n <= 2
End
