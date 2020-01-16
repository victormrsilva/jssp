\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n - 4 x1p + 4 x1n + 4 x2p - 4 x2n + 3 x4p - 3 x4n
Subject To
 constr(0): - 4 x0p + 4 x0n - x1p + x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n
   >= -16
 constr(1): 4 x0p - 4 x0n + 4 x1p - 4 x1n + 4 x2p - 4 x2n - 3 x3p + 3 x3n
   + 3 x4p - 3 x4n >= -8
 constr(2): - 2 x0p + 2 x0n + 3 x1p - 3 x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n
   + x4p - x4n <= 22
 constr(3): - 2 x0p + 2 x0n + 5 x1p - 5 x1n - 3 x2p + 3 x2n - 3 x3p + 3 x3n
   + 3 x4p - 3 x4n = -0
Bounds
 x0p = 0
 x0n <= 4
 x1p <= 2
 x1n <= 4
 x2p <= 3
 x2n <= 4
 x3p <= 1
 x3n <= 3
 x4p <= 3
 x4n <= 3
End
