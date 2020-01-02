\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0 + 2 x1 - 2 x2 - 3 x3
Subject To
 constr(0): 4 x0 + 4 x1 + x2 + 3 x3 <= 12
 constr(1): 3 x0 - 3 x2 - 5 x3 >= -20
 constr(2): - 3 x1 + 2 x2 + 3 x3 = 6
Bounds
 x0 <= 5
 x1 <= 4
 x2 <= 5
 x3 <= 2
Generals
 x0 x1 x2 x3
End
