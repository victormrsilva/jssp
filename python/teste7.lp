\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n + 5 x3p - 5 x3n - 2 x4
Subject To
 constr(0): 4 x1 + 5 x2p - 5 x2n - 3 x3p + 3 x3n + 4 x4 = 2
 constr(1): 5 x0p - 5 x0n - 5 x1 + 2 x2p - 2 x2n - 4 x4 >= -18
 constr(2): - 5 x0p + 5 x0n + 2 x2p - 2 x2n + 3 x3p - 3 x3n - x4 >= -4
 constr(3): - 4 x0p + 4 x0n - 3 x1 - 3 x2p + 3 x2n - 2 x3p + 2 x3n - x4
   >= -22
Bounds
 x0p <= 1
 x0n <= 1
 x1 <= 1
 x2p <= 1
 x2n <= 1
 x3p <= 5
 x3n <= 2
 x4 <= 1
Generals
 x0p x0n x1 x2p x2n x3p x3n x4
End
