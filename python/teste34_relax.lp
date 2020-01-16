\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n + 2 x1p - 2 x1n
Subject To
 constr(0): 5 x1p - 5 x1n <= 3
 constr(1): 3 x0p - 3 x0n + 3 x1p - 3 x1n >= -10
 constr(2): 2 x0p - 2 x0n - 3 x1p + 3 x1n <= 5
 constr(3): - 3 x0p + 3 x0n - 4 x1p + 4 x1n <= 6
Bounds
 x0p <= 2
 x0n <= 3
 x1p = 0
 x1n <= 2
End
