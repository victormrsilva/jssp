\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0 + 5 x1 - 2 x2
Subject To
 constr(0): - 2 x0 + 2 x1 + 3 x2 - 5 x3 >= -0
 constr(1): 2 x0 - 4 x1 + 4 x2 - x3 = 3
 constr(2): - 5 x0 + 5 x2 >= -16
Bounds
 x0 <= 2
 x1 <= 1
 x2 <= 2
 x3 <= 2
Generals
 x0 x1 x2 x3
End
