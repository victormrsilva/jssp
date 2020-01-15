\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x1 + 3 x2p - 3 x2n - 5 x3p + 5 x3n + 5 x4p - 5 x4n
Subject To
 constr(0): 3 x0p - 3 x0n + 3 x1 - 2 x2p + 2 x2n + x3p - x3n + 5 x4p
   - 5 x4n >= -0
 constr(1): 5 x0p - 5 x0n + 4 x1 + 4 x2p - 4 x2n - x3p + x3n - 5 x4p
   + 5 x4n = 13
 constr(2): - 5 x0p + 5 x0n + 3 x1 - x3p + x3n - 5 x4p + 5 x4n >= -0
 constr(3): - 2 x0p + 2 x0n + 3 x1 >= -19
Bounds
 x0p <= 5
 x0n <= 3
 x1 <= 4
 x2p = 0
 x2n <= 1
 x3p <= 3
 x3n <= 3
 x4p <= 3
 x4n <= 1
Generals
 x0p x0n x1 x2p x2n x3p x3n x4p x4n
End
