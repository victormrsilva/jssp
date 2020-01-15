\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n - 2 x1p + 2 x1n - 3 x2
Subject To
 constr(0): 2 x0p - 2 x0n + 3 x1p - 3 x1n - x2 <= -0
 constr(1): 3 x0p - 3 x0n + x1p - x1n + 5 x2 = 6
Bounds
 x0p = 0
 x0n <= 1
 x1p <= 4
 x1n <= 5
 x2 <= 2
Generals
 x0p x0n x1p x1n x2
End
