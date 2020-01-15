\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x2p - 2 x2n + 3 x3p - 3 x3n + 5 x4
Subject To
 constr(0): - 2 x0p + 2 x0n - 5 x1p + 5 x1n + 3 x2p - 3 x2n - 4 x3p + 4 x3n
   + 3 x4 >= -0
 constr(1): 4 x0p - 4 x0n + 3 x1p - 3 x1n + 5 x3p - 5 x3n - 3 x4 >= -0
 constr(2): 4 x1p - 4 x1n + 2 x3p - 2 x3n + 4 x4 = 18
Bounds
 x0p <= 4
 x0n <= 4
 x1p <= 5
 x1n <= 5
 x2p <= 5
 x2n <= 1
 x3p <= 3
 x3n <= 2
 x4 <= 2
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4
End
