\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0p - 3 x0n - 5 x1p + 5 x1n + 2 x3p - 2 x3n
Subject To
 constr(0): - 4 x0p + 4 x0n + x1p - x1n - 3 x2p + 3 x2n - 4 x3p + 4 x3n
   <= 2
 constr(1): - 3 x1p + 3 x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n <= 19
 constr(2): 4 x0p - 4 x0n + 5 x1p - 5 x1n + 4 x2p - 4 x2n - 3 x3p + 3 x3n
   >= -14
 constr(3): - 2 x2p + 2 x2n + 5 x3p - 5 x3n >= -0
 constr(4): - 5 x0p + 5 x0n + 4 x1p - 4 x1n + 3 x2p - 3 x2n - x3p + x3n
   >= -14
Bounds
 x0p <= 1
 x0n <= 1
 x1p = 0
 x1n <= 1
 x2p <= 2
 x2n <= 1
 x3p = 0
 x3n <= 1
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
