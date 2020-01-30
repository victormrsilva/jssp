\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n - 4 x2p + 4 x2n
Subject To
 constr(0): 5 x0p - 5 x0n + 3 x1p - 3 x1n - 3 x2p + 3 x2n - 4 x3p + 4 x3n
   = 2
 constr(1): x0p - x0n + 5 x1p - 5 x1n + 3 x2p - 3 x2n - x3p + x3n >= -0
 constr(2): - 5 x0p + 5 x0n + 2 x1p - 2 x1n + 3 x2p - 3 x2n - 3 x3p + 3 x3n
   >= -0
 constr(3): - 3 x0p + 3 x0n - 3 x1p + 3 x1n + 4 x2p - 4 x2n - x3p + x3n
   = 1
Bounds
 x0p <= 3
 x0n <= 5
 x1p <= 4
 x1n <= 5
 x2p <= 3
 x2n <= 1
 x3p <= 1
 x3n <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
