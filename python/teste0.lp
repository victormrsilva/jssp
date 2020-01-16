\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  4 x0p - 4 x0n + x1 + 3 x2p - 3 x2n + 3 x3p - 3 x3n + 4 x4p - 4 x4n
Subject To
 constr(0): 4 x0p - 4 x0n + x1 - 4 x2p + 4 x2n - 5 x4p + 5 x4n <= 15
 constr(1): - 4 x1 - 4 x2p + 4 x2n + x4p - x4n = 7
 constr(2): 3 x0p - 3 x0n + x2p - x2n + 4 x3p - 4 x3n + 3 x4p - 3 x4n
   <= 17
Bounds
 x0p <= 4
 x0n <= 2
 x1 <= 5
 x2p <= 2
 x2n <= 3
 x3p <= 3
 x3n <= 1
 x4p = 0
 x4n <= 2
Generals
 x0p x0n x1 x2p x2n x3p x3n x4p x4n
End
