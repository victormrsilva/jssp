\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0 + 4 x1p - 4 x1n - 4 x2p + 4 x2n
Subject To
 constr(0): - 2 x0 - x1p + x1n - 4 x2p + 4 x2n + 4 x3p - 4 x3n >= -0
 constr(1): - 4 x0 - 5 x2p + 5 x2n - 5 x3p + 5 x3n = 5
Bounds
 x0 <= 4
 x1p <= 4
 x1n <= 2
 x2p <= 3
 x2n <= 4
 x3p <= 3
 x3n <= 5
Generals
 x0 x1p x1n x2p x2n x3p x3n
End
