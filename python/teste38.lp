\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n + 3 x1p - 3 x1n - 5 x2p + 5 x2n - x3p + x3n + 3 x4p - 3 x4n
Subject To
 constr(0): x0p - x0n - 5 x1p + 5 x1n + 4 x2p - 4 x2n + 4 x3p - 4 x3n + x4p
   - x4n <= 23
 constr(1): - 3 x1p + 3 x1n + 2 x2p - 2 x2n + 4 x3p - 4 x3n <= 4
 constr(2): 4 x0p - 4 x0n + 3 x1p - 3 x1n + 3 x2p - 3 x2n + 4 x4p - 4 x4n
   >= -15
 constr(3): 3 x0p - 3 x0n + x1p - x1n - 2 x3p + 2 x3n >= -0
Bounds
 x0p <= 3
 x0n <= 4
 x1p <= 5
 x1n <= 2
 x2p = 0
 x2n <= 2
 x3p <= 4
 x3n <= 2
 x4p <= 1
 x4n <= 5
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n x4p x4n
End
