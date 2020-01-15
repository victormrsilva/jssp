\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - 3 x1p + 3 x1n + x2p - x2n + x3p - x3n
Subject To
 constr(0): 5 x0p - 5 x0n - 3 x1p + 3 x1n - x2p + x2n - 5 x3p + 5 x3n
   >= -15
 constr(1): - 5 x0p + 5 x0n + 2 x2p - 2 x2n - 5 x3p + 5 x3n <= 15
 constr(2): 4 x1p - 4 x1n - x2p + x2n - 4 x3p + 4 x3n = 4
Bounds
 x0p <= 5
 x0n <= 5
 x1p <= 5
 x1n <= 1
 x2p <= 1
 x2n <= 5
 x3p <= 3
 x3n <= 1
End
