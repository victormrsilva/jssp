\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0p - 3 x0n - 4 x1p + 4 x1n
Subject To
 constr(0): 3 x0p - 3 x0n + 2 x1p - 2 x1n <= 8
 constr(1): - 2 x0p + 2 x0n <= 3
 constr(2): - x0p + x0n - 2 x1p + 2 x1n = -0
 constr(3): x0p - x0n + 2 x1p - 2 x1n <= 9
Bounds
 x0p <= 1
 x0n <= 5
 x1p <= 5
 x1n <= 3
Generals
 x0p x0n x1p x1n
End
