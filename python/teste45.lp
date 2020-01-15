\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - 5 x1p + 5 x1n - 4 x3
Subject To
 constr(0): - x0p + x0n - 5 x1p + 5 x1n + 3 x2p - 3 x2n - 3 x3 <= 4
 constr(1): - 3 x0p + 3 x0n - 5 x2p + 5 x2n + 2 x3 = 5
 constr(2): 5 x0p - 5 x0n + 2 x1p - 2 x1n - 5 x2p + 5 x2n - 5 x3 <= 1
 constr(3): - 5 x0p + 5 x0n + 3 x1p - 3 x1n + x2p - x2n <= 13
Bounds
 x0p = 0
 x0n <= 2
 x1p <= 5
 x1n <= 1
 x2p <= 2
 x2n <= 4
 x3 <= 4
Generals
 x0p x0n x1p x1n x2p x2n x3
End
