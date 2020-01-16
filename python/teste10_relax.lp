\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  3 x1p - 3 x1n + 2 x2p - 2 x2n - 3 x3p + 3 x3n
Subject To
 constr(0): - 3 x0p + 3 x0n + 2 x3p - 2 x3n <= 7
 constr(1): - 5 x0p + 5 x0n + 5 x1p - 5 x1n + 5 x3p - 5 x3n >= -0
 constr(2): - 4 x0p + 4 x0n + x2p - x2n <= 16
 constr(3): - x1p + x1n - 2 x2p + 2 x2n + 4 x3p - 4 x3n <= 19
 constr(4): 4 x0p - 4 x0n - x1p + x1n + 3 x2p - 3 x2n + 3 x3p - 3 x3n
   >= -1
Bounds
 x0p = 0
 x0n <= 3
 x1p <= 3
 x1n <= 5
 x2p <= 2
 x2n <= 5
 x3p <= 5
 x3n <= 4
End
