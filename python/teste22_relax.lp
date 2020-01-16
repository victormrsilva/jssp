\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x1p - x1n - 3 x2p + 3 x2n - 3 x3p + 3 x3n + 3 x4p - 3 x4n
Subject To
 constr(0): 2 x1p - 2 x1n + 4 x3p - 4 x3n + 3 x4p - 3 x4n >= -0
 constr(1): - x0p + x0n + x1p - x1n + 5 x2p - 5 x2n + 2 x3p - 2 x3n + 3 x4p
   - 3 x4n <= 18
Bounds
 x0p = 0
 x0n <= 4
 x1p <= 1
 x1n <= 4
 x2p <= 5
 x2n <= 2
 x3p <= 3
 x3n <= 5
 x4p <= 3
 x4n <= 1
End
