\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0 + 4 x1 - 5 x2p + 5 x2n + 2 x3 + 4 x4p - 4 x4n
Subject To
 constr(0): 2 x1 - x2p + x2n + x3 + 4 x4p - 4 x4n <= 16
 constr(1): 3 x1 + 3 x2p - 3 x2n + 3 x3 - 4 x4p + 4 x4n = 11
 constr(2): - 4 x0 + 2 x2p - 2 x2n + 3 x3 - 5 x4p + 5 x4n >= -0
 constr(3): - 3 x0 - x1 + 5 x2p - 5 x2n - 2 x3 - x4p + x4n >= -0
 constr(4): - 5 x2p + 5 x2n - x3 - 5 x4p + 5 x4n >= -0
Bounds
 x0 <= 1
 x1 <= 5
 x2p <= 2
 x2n <= 2
 x3 <= 4
 x4p = 0
 x4n <= 4
Generals
 x0 x1 x2p x2n x3 x4p x4n
End
