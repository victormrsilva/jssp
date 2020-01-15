\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n + 3 x2p - 3 x2n + 4 x3
Subject To
 constr(0): - x0p + x0n + 3 x1 - 5 x2p + 5 x2n <= 16
 constr(1): 4 x0p - 4 x0n + 3 x1 >= -0
 constr(2): 4 x0p - 4 x0n + 4 x1 - 5 x2p + 5 x2n + x3 <= 8
Bounds
 x0p <= 4
 x0n <= 1
 x1 <= 5
 x2p = 0
 x2n <= 1
 x3 <= 2
End
