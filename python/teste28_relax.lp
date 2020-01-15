\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0p - 3 x0n + 3 x1p - 3 x1n + 4 x2
Subject To
 constr(0): 3 x0p - 3 x0n + 4 x1p - 4 x1n - x2 = 11
 constr(1): - 3 x0p + 3 x0n + x1p - x1n + 2 x2 <= 14
Bounds
 x0p <= 3
 x0n <= 1
 x1p <= 1
 x1n <= 4
 x2 <= 2
End
