\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n - 4 x2p + 4 x2n
Subject To
 constr(0): 5 x0p - 5 x0n - 2 x1p + 2 x1n - x2p + x2n = 8
 constr(1): 3 x0p - 3 x0n + 3 x2p - 3 x2n >= -14
 constr(2): - 5 x0p + 5 x0n - 2 x1p + 2 x1n + 4 x2p - 4 x2n <= 4
Bounds
 x0p <= 5
 x0n <= 2
 x1p <= 3
 x1n <= 2
 x2p <= 4
 x2n <= 5
Generals
 x0p x0n x1p x1n x2p x2n
End
