\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x1 + 4 x2 + x4
Subject To
 constr(0): - 4 x0 + 5 x1 - x2 + 5 x3 <= 6
 constr(1): - 5 x0 - 3 x1 + 5 x2 - 5 x3 - 5 x4 >= -18
 constr(2): 4 x0 + 5 x1 + x2 - x3 - x4 = 17
Bounds
 x0 <= 3
 x1 <= 4
 x2 <= 4
 x3 <= 5
 x4 <= 3
End
