\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  2 x0p - 2 x0n - 2 x1p + 2 x1n
Subject To
 constr(0): x0p - x0n - 4 x2p + 4 x2n = 12
 constr(1): - 3 x0p + 3 x0n + 4 x1p - 4 x1n - 5 x2p + 5 x2n <= 13
 constr(2): x0p - x0n - 2 x2p + 2 x2n >= -4
Bounds
 x0p <= 5
 x0n <= 5
 x1p <= 5
 x1n <= 2
 x2p <= 4
 x2n <= 5
End
