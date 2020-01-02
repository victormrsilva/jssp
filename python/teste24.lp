\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0 + 3 x2
Subject To
 constr(0): 2 x0 - 3 x1 + x2 <= 9
 constr(1): - 2 x0 + 2 x1 = 2
 constr(2): - x0 - 3 x1 >= -7
 constr(3): 2 x0 + x1 + 5 x2 = 4
Bounds
 x0 <= 5
 x1 <= 3
 x2 <= 1
Generals
 x0 x1 x2
End
