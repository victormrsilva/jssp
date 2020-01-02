\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0 + 5 x1 - 2 x2 + 2 x3
Subject To
 constr(0): - 4 x0 - 5 x2 >= -18
 constr(1): 3 x0 - 2 x2 - 3 x3 <= 13
 constr(2): - 4 x0 + 5 x1 + 5 x2 - 3 x3 = 11
 constr(3): x1 + 2 x2 + 4 x3 = 19
Bounds
 x0 <= 1
 x1 <= 5
 x2 <= 5
 x3 <= 3
End
