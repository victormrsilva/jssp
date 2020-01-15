\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n + 3 x1p - 3 x1n - 3 x2p + 3 x2n + 3 x4p - 3 x4n
Subject To
 constr(0): - x0p + x0n + 4 x1p - 4 x1n - 4 x2p + 4 x2n - x4p + x4n >= -0
 constr(1): - 5 x0p + 5 x0n - x1p + x1n + 4 x2p - 4 x2n + 4 x3 + x4p - x4n
   = 12
 constr(2): - 2 x0p + 2 x0n - 3 x1p + 3 x1n + 5 x2p - 5 x2n + 5 x3 - 2 x4p
   + 2 x4n <= 22
Bounds
 x0p <= 3
 x0n <= 2
 x1p = 0
 x1n <= 5
 x2p <= 3
 x2n <= 5
 x3 <= 3
 x4p <= 5
 x4n <= 4
End
