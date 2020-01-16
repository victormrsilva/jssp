\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x1p - 3 x1n + x3
Subject To
 constr(0): 5 x0 + 4 x1p - 4 x1n + 5 x2 + 3 x3 = 20
 constr(1): 2 x0 - 3 x1p + 3 x1n - 5 x2 + 5 x3 <= 10
 constr(2): - 2 x0 + 2 x2 - x3 <= 15
 constr(3): - 5 x0 - 5 x1p + 5 x1n - 4 x2 + 5 x3 >= -16
 constr(4): x0 + 4 x1p - 4 x1n + x3 >= -0
Bounds
 x0 <= 3
 x1p <= 3
 x1n <= 4
 x2 <= 1
 x3 <= 4
Generals
 x0 x1p x1n x2 x3
End
