\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n - 3 x1p + 3 x1n - 5 x2p + 5 x2n
Subject To
 constr(0): - 4 x0p + 4 x0n - 2 x1p + 2 x1n + 3 x2p - 3 x2n <= 10
 constr(1): 3 x0p - 3 x0n + 4 x1p - 4 x1n + 2 x2p - 2 x2n >= -9
 constr(2): - 2 x0p + 2 x0n - 5 x1p + 5 x1n + 3 x2p - 3 x2n >= -0
 constr(3): - 2 x1p + 2 x1n + 5 x2p - 5 x2n <= 2
Bounds
 x0p <= 4
 x0n <= 2
 x1p = 0
 x1n <= 4
 x2p <= 1
 x2n <= 4
Generals
 x0p x0n x1p x1n x2p x2n
End
