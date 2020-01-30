\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - 3 x1
Subject To
 constr(0): - 5 x0p + 5 x0n + 4 x1 + 4 x2p - 4 x2n - 2 x3 <= 15
 constr(1): 5 x1 + 5 x2p - 5 x2n - x3 <= 4
 constr(2): 5 x0p - 5 x0n - 5 x1 + 3 x2p - 3 x2n - x3 <= 5
 constr(3): 2 x1 - 4 x2p + 4 x2n + 4 x3 = 4
Bounds
 x0p <= 2
 x0n <= 4
 x1 <= 1
 x2p <= 1
 x2n <= 1
 x3 <= 3
End
