\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0 - 5 x1 - 2 x3 - 4 x4
Subject To
 constr(0): 3 x0 - 3 x1 + 4 x2 + 3 x3 - 4 x4 <= 19
 constr(1): x2 + 2 x3 + 3 x4 = 3
 constr(2): 2 x0 + 2 x2 + 3 x3 - x4 >= -0
 constr(3): x0 - 4 x1 - 4 x2 - 2 x3 + 3 x4 >= -0
Bounds
 x0 <= 2
 x1 <= 5
 x2 <= 5
 x3 <= 2
 x4 <= 5
Generals
 x0 x1 x2 x3 x4
End
