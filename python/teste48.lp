\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x1 + 3 x2
Subject To
 constr(0): - 2 x0 + 3 x1 + 3 x2 + 4 x3 >= -0
 constr(1): 5 x0 - 3 x1 + 2 x2 + 4 x3 >= -0
 constr(2): x0 + 3 x1 - x2 + 2 x3 <= 13
 constr(3): - 5 x0 + 4 x1 + 3 x2 - 2 x3 <= 19
 constr(4): 3 x1 + x2 + 4 x3 = 14
Bounds
 x0 <= 2
 x1 <= 2
 x2 <= 2
 x3 <= 5
Generals
 x0 x1 x2 x3
End
