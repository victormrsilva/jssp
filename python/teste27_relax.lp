\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - x1 + 2 x2p - 2 x2n - 4 x3p + 4 x3n + 3 x4p - 3 x4n
Subject To
 constr(0): - 4 x0p + 4 x0n + 2 x1 + 5 x3p - 5 x3n <= 19
 constr(1): - 2 x0p + 2 x0n + x1 - 5 x2p + 5 x2n + x3p - x3n = 21
Bounds
 x0p <= 3
 x0n <= 4
 x1 <= 4
 x2p <= 5
 x2n <= 4
 x3p <= 3
 x3n <= 2
 x4p <= 5
 x4n <= 4
End
