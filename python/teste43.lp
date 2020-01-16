\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n
Subject To
 constr(0): x1p - x1n <= 8
 constr(1): x0p - x0n - 2 x1p + 2 x1n <= 5
 constr(2): - x0p + x0n - 4 x1p + 4 x1n = 8
Bounds
 x0p <= 1
 x0n <= 4
 x1p <= 4
 x1n <= 3
Generals
 x0p x0n x1p x1n
End
