\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n - 4 x1p + 4 x1n
Subject To
 constr(0): 3 x0p - 3 x0n + 4 x1p - 4 x1n <= 9
 constr(1): - x0p + x0n - 2 x1p + 2 x1n <= 10
Bounds
 x0p <= 2
 x0n <= 3
 x1p <= 1
 x1n <= 2
End
