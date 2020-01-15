\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0 + 2 x1 - 3 x2p + 3 x2n + x3p - x3n + 4 x4p - 4 x4n
Subject To
 constr(0): 4 x0 - 3 x1 - x2p + x2n - 2 x3p + 2 x3n - 4 x4p + 4 x4n <= 14
 constr(1): - 4 x1 + 2 x2p - 2 x2n - x3p + x3n + 2 x4p - 2 x4n = 4
 constr(2): - 4 x0 + x1 - 3 x2p + 3 x2n - 4 x3p + 4 x3n - 2 x4p + 2 x4n
   = 16
 constr(3): - 4 x0 - 5 x1 - 5 x2p + 5 x2n + 3 x3p - 3 x3n + x4p - x4n <= 4
Bounds
 x0 <= 2
 x1 <= 4
 x2p = 0
 x2n <= 4
 x3p <= 5
 x3n <= 5
 x4p = 0
 x4n <= 5
End
