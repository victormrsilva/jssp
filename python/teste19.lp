\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n + 2 x1p - 2 x1n + 4 x2p - 4 x2n - x3p + x3n
Subject To
 constr(0): - 3 x0p + 3 x0n - 2 x1p + 2 x1n - 5 x2p + 5 x2n - 5 x3p + 5 x3n
   <= 19
 constr(1): 2 x1p - 2 x1n + 4 x2p - 4 x2n + 3 x3p - 3 x3n = 10
 constr(2): - 5 x0p + 5 x0n + 4 x1p - 4 x1n - x3p + x3n >= -2
Bounds
 x0p <= 1
 x0n <= 3
 x1p = 0
 x1n <= 2
 x2p <= 5
 x2n <= 1
 x3p <= 3
 x3n <= 3
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
