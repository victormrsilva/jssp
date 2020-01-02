\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0 - x1 - 4 x2 + 3 x3 - 3 x4
Subject To
 constr(0): - 3 x0 + 5 x1 + 5 x3 + 3 x4 <= 14
 constr(1): - 3 x0 - 5 x2 - 2 x3 + 5 x4 <= 25
 constr(2): - 5 x0 + 2 x1 - 5 x2 - x4 = 2
 constr(3): - 3 x0 + 4 x1 + 2 x2 + 4 x3 - 3 x4 >= -0
Bounds
 x0 <= 2
 x1 <= 4
 x2 <= 3
 x3 <= 5
 x4 <= 3
End
