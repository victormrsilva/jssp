\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n
Subject To
 constr(0): - 3 x0p + 3 x0n + 3 x1 <= 4
 constr(1): 4 x0p - 4 x0n + 5 x1 <= 3
 constr(2): 2 x0p - 2 x0n - x1 <= 7
Bounds
 x0p <= 5
 x0n <= 3
 x1 <= 5
Generals
 x0p x0n x1
End
