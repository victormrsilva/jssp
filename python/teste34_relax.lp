\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0 - 4 x1 - 2 x2 - 3 x3 + 2 x4
Subject To
 constr(0): - x0 + 4 x1 - 4 x2 - 3 x3 - 2 x4 >= -1
 constr(1): - 5 x1 - 3 x2 - 3 x3 - 3 x4 >= -12
 constr(2): 3 x0 + 5 x1 + x2 - x3 - 5 x4 <= 6
 constr(3): 2 x0 - 5 x1 - 2 x2 - 3 x3 - 3 x4 <= 7
 constr(4): 3 x0 - 5 x1 + x2 + 5 x3 + 5 x4 >= -0
Bounds
 x0 <= 2
 x1 <= 2
 x2 <= 2
 x3 <= 3
 x4 <= 3
End
