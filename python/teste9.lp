\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n - 5 x1p + 5 x1n - 3 x2p + 3 x2n + x3 - 3 x4p + 3 x4n
Subject To
 constr(0): 3 x0p - 3 x0n - 3 x1p + 3 x1n - x2p + x2n - 2 x3 - 5 x4p
   + 5 x4n >= -7
 constr(1): - 2 x0p + 2 x0n - 4 x1p + 4 x1n + 3 x2p - 3 x2n - 5 x3 - 5 x4p
   + 5 x4n <= 5
Bounds
 x0p <= 2
 x0n <= 5
 x1p <= 2
 x1n <= 2
 x2p = 0
 x2n <= 1
 x3 <= 3
 x4p <= 2
 x4n <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3 x4p x4n
End
