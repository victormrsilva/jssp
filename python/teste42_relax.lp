\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0 + 5 x1p - 5 x1n
Subject To
 constr(0): 3 x1p - 3 x1n <= 1
 constr(1): - 3 x1p + 3 x1n <= 8
Bounds
 x0 <= 5
 x1p <= 5
 x1n <= 4
End
