\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n + 4 x1p - 4 x1n
Subject To
 constr(0): - 4 x0p + 4 x0n - 3 x1p + 3 x1n <= 7
 constr(1): - 5 x0p + 5 x0n + x1p - x1n = -0
Bounds
 x0p <= 1
 x0n <= 1
 x1p <= 4
 x1n <= 3
Generals
 x0p x0n x1p x1n
End
