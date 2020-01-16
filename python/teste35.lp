\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 2 x0p + 2 x0n + 4 x2p - 4 x2n + 5 x3
Subject To
 constr(0): - 4 x1 - x2p + x2n + 5 x3 = 10
 constr(1): 3 x0p - 3 x0n + 5 x1 - 4 x2p + 4 x2n + 3 x3 = -0
 constr(2): - 2 x0p + 2 x0n - 2 x1 + 3 x2p - 3 x2n - 5 x3 <= 10
Bounds
 x0p <= 4
 x0n <= 4
 x1 <= 2
 x2p <= 2
 x2n <= 2
 x3 <= 4
Generals
 x0p x0n x1 x2p x2n x3
End
