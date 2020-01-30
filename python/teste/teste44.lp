\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0p - x0n + 3 x2p - 3 x2n
Subject To
 constr(0): - 4 x0p + 4 x0n + 2 x1p - 2 x1n = 8
 constr(1): - x1p + x1n - x2p + x2n >= -1
 constr(2): - 3 x0p + 3 x0n - x1p + x1n - 5 x2p + 5 x2n <= 14
Bounds
 x0p = 0
 x0n <= 5
 x1p <= 3
 x1n <= 3
 x2p <= 4
 x2n <= 1
Generals
 x0p x0n x1p x1n x2p x2n
End
