\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0 - x1p + x1n
Subject To
 constr(0): 2 x0 - 5 x1p + 5 x1n + 2 x3p - 2 x3n >= -0
 constr(1): - 4 x0 + 5 x1p - 5 x1n + 2 x2p - 2 x2n - 2 x3p + 2 x3n <= 10
 constr(2): - 3 x0 - 2 x1p + 2 x1n - 3 x2p + 3 x2n <= 11
 constr(3): 2 x0 + x1p - x1n + 3 x2p - 3 x2n + 3 x3p - 3 x3n = 8
Bounds
 x0 <= 4
 x1p <= 5
 x1n <= 2
 x2p <= 5
 x2n <= 3
 x3p = 0
 x3n <= 5
End
