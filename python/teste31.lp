\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n - 2 x2p + 2 x2n
Subject To
 constr(0): - 4 x0p + 4 x0n + 3 x1p - 3 x1n + 4 x2p - 4 x2n <= 4
 constr(1): - 2 x0p + 2 x0n + x1p - x1n - 2 x2p + 2 x2n >= -9
Bounds
 x0p <= 4
 x0n <= 1
 x1p = 0
 x1n <= 5
 x2p <= 1
 x2n <= 1
Generals
 x0p x0n x1p x1n x2p x2n
End
