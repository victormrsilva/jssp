\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0p + 4 x0n + x1p - x1n
Subject To
 constr(0): 5 x0p - 5 x0n + 2 x1p - 2 x1n <= 2
 constr(1): - 3 x0p + 3 x0n - 5 x1p + 5 x1n = 10
Bounds
 x0p <= 1
 x0n <= 2
 x1p <= 2
 x1n <= 5
End
