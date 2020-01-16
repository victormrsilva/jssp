\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x1p + 4 x1n
Subject To
 constr(0): 2 x0 + 4 x1p - 4 x1n <= 3
 constr(1): - 3 x0 + 2 x1p - 2 x1n >= -3
Bounds
 x0 <= 3
 x1p <= 2
 x1n <= 1
End
