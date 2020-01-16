\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0p + 4 x0n + 5 x1p - 5 x1n - 2 x2p + 2 x2n
Subject To
 constr(0): 3 x0p - 3 x0n - 4 x2p + 4 x2n <= 4
 constr(1): 2 x1p - 2 x1n = -0
 constr(2): - 4 x0p + 4 x0n + 3 x1p - 3 x1n - 3 x2p + 3 x2n <= 12
 constr(3): - 2 x0p + 2 x0n - x2p + x2n <= 6
 constr(4): - 3 x0p + 3 x0n - 3 x1p + 3 x1n - 5 x2p + 5 x2n >= -0
Bounds
 x0p <= 3
 x0n <= 5
 x1p <= 3
 x1n <= 3
 x2p <= 2
 x2n <= 3
Generals
 x0p x0n x1p x1n x2p x2n
End
