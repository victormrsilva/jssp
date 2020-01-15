\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n - 4 x3p + 4 x3n + 2 x4p - 2 x4n
Subject To
 constr(0): - x0p + x0n + 3 x1p - 3 x1n - 4 x2p + 4 x2n - 2 x3p + 2 x3n
   - 5 x4p + 5 x4n <= 2
 constr(1): - 4 x0p + 4 x0n + 5 x1p - 5 x1n + 5 x2p - 5 x2n - 4 x4p + 4 x4n
   >= -16
 constr(2): 5 x0p - 5 x0n - 4 x1p + 4 x1n + 2 x2p - 2 x2n - 2 x3p + 2 x3n
   - 2 x4p + 2 x4n = 7
Bounds
 x0p = 0
 x0n <= 1
 x1p <= 5
 x1n <= 1
 x2p <= 5
 x2n <= 1
 x3p = 0
 x3n <= 2
 x4p <= 4
 x4n <= 5
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4p x4n
End
