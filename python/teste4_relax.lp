\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n + 5 x1p - 5 x1n
Subject To
 constr(0): - 5 x0p + 5 x0n + x1p - x1n = 5
 constr(1): - 4 x0p + 4 x0n <= 9
 constr(2): - x0p + x0n - 5 x1p + 5 x1n <= 9
Bounds
 x0p = 0
 x0n <= 3
 x1p <= 2
 x1n <= 1
End
