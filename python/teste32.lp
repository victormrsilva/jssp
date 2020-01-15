\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n + x1p - x1n - 2 x2p + 2 x2n - 4 x3
Subject To
 constr(0): 5 x0p - 5 x0n + 4 x1p - 4 x1n - 3 x2p + 3 x2n = 12
 constr(1): - 5 x0p + 5 x0n + 5 x1p - 5 x1n + 4 x2p - 4 x2n + 5 x3 = 13
 constr(2): - 4 x0p + 4 x0n + 3 x1p - 3 x1n + 5 x2p - 5 x2n + 2 x3 <= 11
 constr(3): - 2 x0p + 2 x0n - 4 x1p + 4 x1n - 2 x2p + 2 x2n + x3 >= -17
 constr(4): - 3 x0p + 3 x0n - 5 x1p + 5 x1n - 2 x2p + 2 x2n - x3 >= -10
Bounds
 x0p <= 1
 x0n <= 2
 x1p <= 5
 x1n <= 4
 x2p <= 3
 x2n <= 4
 x3 <= 2
Generals
 x0p x0n x1p x1n x2p x2n x3
End
