\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - x0p + x0n + 2 x1p - 2 x1n - 4 x2p + 4 x2n
Subject To
 constr(0): - 2 x0p + 2 x0n + 3 x1p - 3 x1n + x2p - x2n >= -0
 constr(1): - 5 x0p + 5 x0n - 3 x1p + 3 x1n + 5 x2p - 5 x2n <= 14
 constr(2): 2 x0p - 2 x0n + 3 x1p - 3 x1n <= 6
 constr(3): - 2 x0p + 2 x0n + 4 x2p - 4 x2n = 2
Bounds
 x0p <= 2
 x0n <= 3
 x1p <= 2
 x1n <= 3
 x2p <= 5
 x2n <= 2
Generals
 x0p x0n x1p x1n x2p x2n
End
