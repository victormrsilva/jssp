\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x2p - 5 x2n - 5 x3
Subject To
 constr(0): 5 x0p - 5 x0n - 3 x1p + 3 x1n + 4 x2p - 4 x2n + 2 x3 = 3
 constr(1): 3 x0p - 3 x0n - x1p + x1n + 3 x2p - 3 x2n <= 9
 constr(2): - x0p + x0n + x1p - x1n + 2 x2p - 2 x2n >= -0
Bounds
 x0p <= 1
 x0n <= 2
 x1p <= 1
 x1n <= 4
 x2p <= 5
 x2n <= 1
 x3 <= 4
End
