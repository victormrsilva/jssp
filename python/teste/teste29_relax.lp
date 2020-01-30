\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n + x1p - x1n
Subject To
 constr(0): - 3 x0p + 3 x0n - 4 x1p + 4 x1n <= 6
 constr(1): - 4 x1p + 4 x1n <= 7
 constr(2): 3 x0p - 3 x0n - 3 x1p + 3 x1n >= -1
Bounds
 x0p <= 2
 x0n <= 3
 x1p <= 3
 x1n <= 2
End
