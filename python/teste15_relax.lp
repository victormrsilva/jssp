\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0 + x1p - x1n - 4 x3p + 4 x3n
Subject To
 constr(0): 3 x1p - 3 x1n + 3 x2p - 3 x2n - 5 x3p + 5 x3n = 19
 constr(1): - x0 - x1p + x1n - 2 x3p + 2 x3n <= 19
Bounds
 x0 <= 5
 x1p <= 2
 x1n <= 4
 x2p <= 3
 x2n <= 4
 x3p = 0
 x3n <= 5
End
