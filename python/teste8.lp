\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0 + 2 x2 + 3 x3
Subject To
 constr(0): 3 x0 - 4 x2 - 2 x3 >= -4
 constr(1): x0 - x1 + 3 x3 >= -0
 constr(2): 5 x0 + 4 x1 - 2 x2 = 5
Bounds
 x0 <= 2
 x1 <= 2
 x2 <= 3
 x3 <= 5
Generals
 x0 x1 x2 x3
End
