\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0 + 5 x1 - 3 x2
Subject To
 constr(0): 2 x0 + 3 x2 = 2
 constr(1): 5 x0 + 5 x1 + x2 + 5 x3 <= 7
 constr(2): x0 - 2 x1 - x2 + 4 x3 <= 8
 constr(3): - 3 x0 - 2 x2 - 3 x3 <= 17
 constr(4): 5 x0 + 4 x1 - 4 x2 >= -14
Bounds
 x0 <= 5
 x1 <= 2
 x2 <= 3
 x3 <= 3
End
