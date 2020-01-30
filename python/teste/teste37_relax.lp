\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - 3 x1p + 3 x1n - 4 x2
Subject To
 constr(0): - 3 x0p + 3 x0n + x1p - x1n + 3 x2 <= 12
 constr(1): 2 x0p - 2 x0n - 4 x1p + 4 x1n + x2 = -0
 constr(2): - 3 x0p + 3 x0n + 2 x2 <= 1
Bounds
 x0p <= 5
 x0n <= 4
 x1p <= 3
 x1n <= 2
 x2 <= 5
End
