\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n - 2 x1p + 2 x1n + x3p - x3n - 3 x4p + 3 x4n
Subject To
 constr(0): 5 x0p - 5 x0n - x2p + x2n + x3p - x3n >= -9
 constr(1): - x0p + x0n + 2 x1p - 2 x1n - 2 x2p + 2 x2n + 5 x3p - 5 x3n
   - 2 x4p + 2 x4n = 25
Bounds
 x0p <= 3
 x0n <= 3
 x1p = 0
 x1n <= 2
 x2p <= 5
 x2n <= 2
 x3p <= 3
 x3n <= 4
 x4p <= 1
 x4n <= 2
End
