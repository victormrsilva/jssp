\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n + 3 x1p - 3 x1n + 3 x2p - 3 x2n - 3 x3 - 2 x4p + 2 x4n
Subject To
 constr(0): 2 x1p - 2 x1n - 5 x2p + 5 x2n - x3 >= -0
 constr(1): 3 x0p - 3 x0n + x1p - x1n + 3 x3 - x4p + x4n = 9
 constr(2): - x1p + x1n - 4 x2p + 4 x2n + 4 x3 + 4 x4p - 4 x4n = 6
Bounds
 x0p <= 2
 x0n <= 1
 x1p <= 4
 x1n <= 1
 x2p <= 2
 x2n <= 2
 x3 <= 4
 x4p <= 4
 x4n <= 5
End
