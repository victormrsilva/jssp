\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0 + 4 x1 - x2 - 3 x3
Subject To
 constr(0): - 2 x0 - 3 x1 + 5 x2 + 5 x3 <= 7
 constr(1): - 2 x2 - 2 x3 <= 6
 constr(2): - x0 + x1 - 2 x2 + 2 x3 >= -3
 constr(3): x0 - 3 x1 + 5 x2 + 3 x3 = 3
 constr(4): 4 x0 + 5 x1 = 14
Bounds
 x0 <= 5
 x1 <= 5
 x2 <= 3
 x3 <= 1
Generals
 x0 x1 x2 x3
End
