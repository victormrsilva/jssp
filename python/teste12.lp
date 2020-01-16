\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n + 5 x1p - 5 x1n - 5 x2p + 5 x2n - 2 x3p + 2 x3n + x4
Subject To
 constr(0): 3 x0p - 3 x0n - 5 x2p + 5 x2n - x3p + x3n + 4 x4 = 4
 constr(1): 3 x1p - 3 x1n + x2p - x2n + 4 x3p - 4 x3n - 2 x4 >= -0
 constr(2): x0p - x0n + 4 x1p - 4 x1n - 5 x2p + 5 x2n - 5 x3p + 5 x3n
   + 4 x4 <= 8
Bounds
 x0p <= 5
 x0n <= 1
 x1p <= 3
 x1n <= 4
 x2p <= 5
 x2n <= 1
 x3p = 0
 x3n <= 1
 x4 <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4
End
