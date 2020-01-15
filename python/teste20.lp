\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0p + 4 x0n - 2 x1p + 2 x1n + x2p - x2n
Subject To
 constr(0): - 2 x1p + 2 x1n >= -11
 constr(1): 3 x0p - 3 x0n - 2 x1p + 2 x1n - 4 x2p + 4 x2n = 9
 constr(2): x0p - x0n - 3 x1p + 3 x1n - 2 x2p + 2 x2n >= -0
 constr(3): 4 x1p - 4 x1n + 5 x2p - 5 x2n <= 12
 constr(4): - 4 x0p + 4 x0n + 3 x2p - 3 x2n <= -0
Bounds
 x0p <= 4
 x0n <= 5
 x1p <= 4
 x1n <= 2
 x2p = 0
 x2n <= 5
Generals
 x0p x0n x1p x1n x2p x2n
End
