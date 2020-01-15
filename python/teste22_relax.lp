\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n - 4 x1p + 4 x1n + 3 x4p - 3 x4n
Subject To
 constr(0): - x0p + x0n + 4 x1p - 4 x1n + 2 x2p - 2 x2n + 4 x3p - 4 x3n
   + 5 x4p - 5 x4n >= -0
 constr(1): 4 x0p - 4 x0n - 3 x1p + 3 x1n + x2p - x2n - 2 x3p + 2 x3n
   + 4 x4p - 4 x4n = 20
 constr(2): - 4 x0p + 4 x0n + 2 x1p - 2 x1n + x2p - x2n + x3p - x3n - x4p
   + x4n >= -0
 constr(3): 3 x0p - 3 x0n - 5 x2p + 5 x2n + x3p - x3n + 3 x4p - 3 x4n
   <= -0
 constr(4): 2 x0p - 2 x0n - 5 x1p + 5 x1n + x2p - x2n + 5 x3p - 5 x3n
   - 4 x4p + 4 x4n >= -0
Bounds
 x0p <= 1
 x0n <= 4
 x1p <= 4
 x1n <= 5
 x2p <= 5
 x2n <= 3
 x3p <= 1
 x3n <= 1
 x4p <= 4
 x4n <= 2
End
