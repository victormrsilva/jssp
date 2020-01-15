\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n + x1p - x1n
Subject To
 constr(0): 3 x0p - 3 x0n + 5 x1p - 5 x1n >= -0
 constr(1): 5 x0p - 5 x0n + 4 x1p - 4 x1n = 2
Bounds
 x0p <= 3
 x0n <= 2
 x1p <= 3
 x1n <= 1
End
