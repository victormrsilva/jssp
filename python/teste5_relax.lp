\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + 4 x2p - 4 x2n
Subject To
 constr(0): - 2 x0p + 2 x0n - 5 x2p + 5 x2n >= -13
 constr(1): 2 x0p - 2 x0n + 5 x1p - 5 x1n + 2 x2p - 2 x2n = 6
 constr(2): - 5 x0p + 5 x0n - 3 x1p + 3 x1n - x2p + x2n <= 5
 constr(3): 4 x0p - 4 x0n + 3 x1p - 3 x1n - 2 x2p + 2 x2n >= -0
Bounds
 x0p = 0
 x0n <= 3
 x1p <= 5
 x1n <= 2
 x2p <= 3
 x2n <= 5
End
