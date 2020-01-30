\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n + 2 x1 + x2p - x2n
Subject To
 constr(0): - 3 x0p + 3 x0n + 5 x1 <= 13
 constr(1): x0p - x0n - 2 x1 - x2p + x2n >= -0
 constr(2): - 2 x0p + 2 x0n - 5 x1 + 3 x2p - 3 x2n >= -4
 constr(3): - 2 x0p + 2 x0n + 2 x2p - 2 x2n <= 14
Bounds
 x0p = 0
 x0n <= 3
 x1 <= 3
 x2p = 0
 x2n <= 5
End
