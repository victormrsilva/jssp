\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - 4 x2p + 4 x2n + 3 x3p - 3 x3n + 5 x4
Subject To
 constr(0): - 5 x0p + 5 x0n + 4 x1p - 4 x1n + 4 x2p - 4 x2n + 2 x3p - 2 x3n
   + 2 x4 = 21
 constr(1): - 4 x0p + 4 x0n + 2 x2p - 2 x2n + 5 x3p - 5 x3n + 5 x4 >= -0
 constr(2): 4 x0p - 4 x0n + 4 x1p - 4 x1n - 3 x2p + 3 x2n + x3p - x3n
   + 3 x4 <= 3
 constr(3): 4 x0p - 4 x0n + 2 x1p - 2 x1n + 5 x2p - 5 x2n + 4 x3p - 4 x3n
   + x4 >= -0
Bounds
 x0p <= 1
 x0n <= 5
 x1p <= 1
 x1n <= 4
 x2p <= 1
 x2n <= 5
 x3p <= 5
 x3n <= 1
 x4 <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4
End
