\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + x2 - 3 x3
Subject To
 constr(0): - 5 x0p + 5 x0n + 2 x1p - 2 x1n + x2 + 2 x3 >= -0
 constr(1): - 3 x1p + 3 x1n + x2 - 4 x3 = 2
 constr(2): - x0p + x0n - 5 x1p + 5 x1n - 4 x2 - x3 >= -16
 constr(3): x2 + 2 x3 = 11
 constr(4): - 5 x0p + 5 x0n + x3 <= 8
Bounds
 x0p <= 1
 x0n <= 5
 x1p <= 4
 x1n <= 5
 x2 <= 3
 x3 <= 4
End
