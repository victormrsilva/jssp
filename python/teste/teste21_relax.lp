\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x1 + x2 - 2 x4
Subject To
 constr(0): - x0 + 4 x1 + 5 x2 + 2 x3 + 4 x4 >= -17
 constr(1): 5 x0 - 3 x1 + 4 x2 - x3 = 2
 constr(2): - 2 x0 - 2 x2 + x3 + 4 x4 <= 5
Bounds
 x0 <= 5
 x1 <= 2
 x2 <= 4
 x3 <= 4
 x4 <= 5
End
