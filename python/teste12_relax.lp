\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - x1p + x1n - x2p + x2n
Subject To
 constr(0): - x0p + x0n - x2p + x2n = 2
 constr(1): - 2 x0p + 2 x0n + 2 x1p - 2 x1n + 5 x2p - 5 x2n <= -0
Bounds
 x0p <= 4
 x0n <= 5
 x1p <= 5
 x1n <= 1
 x2p <= 4
 x2n <= 4
End
