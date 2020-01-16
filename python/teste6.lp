\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n - 4 x1 + x2p - x2n
Subject To
 constr(0): 5 x0p - 5 x0n + 5 x1 + x2p - x2n = 15
 constr(1): - 2 x0p + 2 x0n - 3 x1 + 3 x2p - 3 x2n >= -11
 constr(2): 2 x0p - 2 x0n + x1 + 3 x2p - 3 x2n >= -11
Bounds
 x0p <= 4
 x0n <= 4
 x1 <= 3
 x2p <= 2
 x2n <= 3
Generals
 x0p x0n x1 x2p x2n
End
