\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n + 5 x1 + 5 x2
Subject To
 constr(0): x1 + 5 x2 >= -0
 constr(1): - 3 x0p + 3 x0n - 3 x1 + 4 x2 = 4
 constr(2): 4 x0p - 4 x0n + x1 - 2 x2 <= 10
 constr(3): - 2 x0p + 2 x0n - 2 x1 - 5 x2 <= -0
Bounds
 x0p <= 2
 x0n <= 2
 x1 <= 3
 x2 <= 3
End
