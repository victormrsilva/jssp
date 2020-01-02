\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0 + 3 x1 + 3 x2
Subject To
 constr(0): - 3 x0 - 2 x1 - 3 x2 >= -10
 constr(1): 2 x0 + 4 x1 - 2 x2 = 14
 constr(2): 5 x0 - 2 x1 <= 11
 constr(3): - 2 x0 - 5 x1 - x2 <= 4
 constr(4): 3 x0 - 4 x1 >= -14
Bounds
 x0 <= 5
 x1 <= 3
 x2 <= 5
Generals
 x0 x1 x2
End
