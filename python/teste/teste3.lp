\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0 + 5 x1 + 3 x2
Subject To
 constr(0): 5 x0 + 2 x1 + 2 x2 = 7
 constr(1): - 2 x0 + 3 x1 <= 6
 constr(2): 3 x1 + 3 x2 <= 6
 constr(3): 5 x0 + 3 x2 >= -5
Bounds
 x0 <= 3
 x1 <= 3
 x2 <= 1
Generals
 x0 x1 x2
End
