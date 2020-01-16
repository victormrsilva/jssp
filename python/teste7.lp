\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + x1p - x1n - 4 x2p + 4 x2n + 2 x3p - 2 x3n
Subject To
 constr(0): - x0p + x0n + x1p - x1n - 5 x2p + 5 x2n - x3p + x3n >= -0
 constr(1): - 2 x0p + 2 x0n - x1p + x1n - 5 x2p + 5 x2n + 4 x3p - 4 x3n
   = 1
 constr(2): - 4 x0p + 4 x0n - 2 x1p + 2 x1n + 5 x2p - 5 x2n - 2 x3p + 2 x3n
   >= -20
 constr(3): 3 x0p - 3 x0n + 4 x1p - 4 x1n - 5 x2p + 5 x2n + 3 x3p - 3 x3n
   >= -7
Bounds
 x0p <= 3
 x0n <= 3
 x1p <= 4
 x1n <= 1
 x2p <= 5
 x2n <= 2
 x3p <= 1
 x3n <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
