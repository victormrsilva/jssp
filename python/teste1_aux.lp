\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x0p - 3 x0n - 5 x1p + 5 x1n + 4 x3p - 4 x3n
Subject To
 constr(0): - 5 x0p + 5 x0n + 3 x2p - 3 x2n - 2 x3p + 2 x3n <= 17
 constr(1): - 4 x0p + 4 x0n - x1p + x1n + 3 x2p - 3 x2n = 19
Bounds
 x0p <= 5
 x0n <= 3
 x1p <= 1
 x1n <= 2
 x2p <= 5
 x2n <= 2
 x3p <= 2
 x3n <= 5
Generals
 x0p x0n x1p x1n x2p x2n x3p x3n
End
