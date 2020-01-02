\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0 - 5 x1 - x2 - 5 x3
Subject To
 constr(0): x0 - 3 x1 + 2 x3 <= 15
 constr(1): - 3 x1 - 5 x3 >= -5
 constr(2): - 3 x0 - 3 x1 - 2 x2 - 5 x3 <= 4
Bounds
 x0 <= 2
 x1 <= 2
 x2 <= 5
 x3 <= 4
End
