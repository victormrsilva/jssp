\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0 + x1 - 5 x2 + 2 x3 - 3 x4
Subject To
 constr(0): 5 x0 + 5 x1 - 3 x2 - 5 x3 - 5 x4 <= 14
 constr(1): 2 x1 + 3 x2 + x3 + 4 x4 = 8
Bounds
 -4 <= x0 <= 4
 -3 <= x1 <= 4
 -3 <= x2 <= 5
 x3 <= 3
 -1 <= x4 <= 5
Generals
 x0 x1 x2 x3 x4
End
