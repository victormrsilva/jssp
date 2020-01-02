\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0 + 5 x1 + 2 x2
Subject To
 constr(0): 4 x0 + x1 + 4 x2 + 5 x3 + 2 x4 = 1
 constr(1): - 5 x0 + 5 x1 + x2 + 5 x3 - x4 >= -0
Bounds
 x0 <= 5
 x1 <= 1
 x2 <= 1
 x3 <= 4
 x4 <= 2
Generals
 x0 x1 x2 x3 x4
End
