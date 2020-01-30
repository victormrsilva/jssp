\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x1p - x1n
Subject To
 constr(0): - x0p + x0n - 2 x1p + 2 x1n <= 2
 constr(1): 4 x1p - 4 x1n >= -2
 constr(2): - 2 x0p + 2 x0n + x1p - x1n = 5
 constr(3): - 5 x1p + 5 x1n >= -8
Bounds
 x0p <= 2
 x0n <= 3
 x1p <= 3
 x1n <= 4
Generals
 x0p x0n x1p x1n
End
