\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0p + 4 x0n - 5 x3p + 5 x3n + 3 x4p - 3 x4n
Subject To
 constr(0): - 2 x0p + 2 x0n - x1 + 2 x2p - 2 x2n + 3 x3p - 3 x3n - 2 x4p
   + 2 x4n = -0
 constr(1): - x0p + x0n + 5 x1 - x2p + x2n + 2 x3p - 2 x3n + 4 x4p - 4 x4n
   <= -0
 constr(2): 4 x0p - 4 x0n - x1 - 5 x2p + 5 x2n + 2 x3p - 2 x3n - x4p + x4n
   <= 14
 constr(3): - x0p + x0n + 3 x1 - x2p + x2n + 2 x4p - 2 x4n <= 23
Bounds
 x0p <= 1
 x0n <= 5
 x1 <= 1
 x2p <= 5
 x2n <= 3
 x3p <= 2
 x3n <= 5
 x4p <= 2
 x4n <= 5
End
