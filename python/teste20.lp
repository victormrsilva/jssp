\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0 + 2 x1 - 3 x3
Subject To
 constr(0): 4 x0 + 4 x1 + 4 x2 + 5 x3 <= 6
 constr(1): 3 x1 + 5 x2 - 5 x3 <= 15
 constr(2): - 2 x0 - 5 x1 <= 7
 constr(3): 4 x2 - x3 <= 10
 constr(4): x1 - x2 + 3 x3 <= 19
Bounds
 x0 <= 3
 x1 <= 2
 x2 <= 3
 x3 <= 5
Generals
 x0 x1 x2 x3
End
