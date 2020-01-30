\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - x2p + x2n + x3p - x3n
Subject To
 constr(0): x0p - x0n - 4 x1p + 4 x1n - x2p + x2n + 5 x3p - 5 x3n >= -0
 constr(1): 3 x0p - 3 x0n + 3 x1p - 3 x1n + x2p - x2n + 4 x4 <= 18
 constr(2): - 3 x1p + 3 x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n + 5 x4 <= 3
 constr(3): x0p - x0n + 5 x1p - 5 x1n + 5 x2p - 5 x2n - 5 x3p + 5 x3n
   + 4 x4 <= 13
Bounds
 x0p <= 2
 x0n <= 2
 x1p <= 4
 x1n <= 1
 x2p = 0
 x2n <= 2
 x3p <= 1
 x3n <= 4
 x4 <= 4
End
