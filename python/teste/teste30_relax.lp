\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + 4 x1 - 5 x2p + 5 x2n
Subject To
 constr(0): - 2 x1 + 2 x2p - 2 x2n <= 7
 constr(1): 5 x1 - 3 x2p + 3 x2n <= 1
 constr(2): - x0p + x0n - 3 x1 <= 14
 constr(3): 2 x0p - 2 x0n - 3 x1 + 2 x2p - 2 x2n <= 10
 constr(4): - 3 x0p + 3 x0n + 3 x1 - 3 x2p + 3 x2n >= -0
Bounds
 x0p <= 5
 x0n <= 4
 x1 <= 5
 x2p <= 4
 x2n <= 2
End
