\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0 + 2 x1 + 3 x3 - 5 x4
Subject To
 constr(0): 5 x0 - x1 + 2 x3 + x4 <= 19
 constr(1): 2 x1 - 5 x3 - 2 x4 >= -14
 constr(2): x0 + 3 x2 + x3 - 2 x4 = 5
 constr(3): 5 x0 - 5 x1 - 4 x2 - x3 - 3 x4 >= -0
Bounds
 x0 <= 5
 x1 <= 1
 x2 <= 3
 x3 <= 2
 x4 <= 3
End
