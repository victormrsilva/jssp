\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0 + 5 x1 - x2 + x3
Subject To
 constr(0): - 4 x0 + 4 x1 + 4 x2 - 3 x3 = 9
 constr(1): - 2 x1 - 2 x3 >= -20
Bounds
 x0 <= 3
 x1 <= 2
 x2 <= 5
 x3 <= 5
End
