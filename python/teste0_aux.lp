\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n + x1p - x1n - 5 x2p + 5 x2n + 2 x3 - 3 x4p + 3 x4n
Subject To
 constr(0): 5 x0p - 5 x0n + 5 x1p - 5 x1n - 3 x2p + 3 x2n - 5 x3 - 5 x4p + 5 x4n <= 14
 constr(1): 2 x1p - 2 x1n + 3 x2p - 3 x2n + x3 + 4 x4p - 4 x4n = 8
Bounds
 x0p <= 4
 x0n <= 4
 x1p <= 4
 x1n <= 3
 x2p <= 5
 x2n <= 3
 x3 <= 3
 x4p <= 5
 x4n <= 1
Generals
 x0p x0n x1p x1n x2p x2n x3 x4p x4n
End
