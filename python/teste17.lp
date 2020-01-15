\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - 4 x1 + 5 x2p - 5 x2n + 4 x3
Subject To
 constr(0): - 4 x0p + 4 x0n + 5 x1 - 2 x2p + 2 x2n - 3 x3 - 3 x4p + 3 x4n
   >= -0
 constr(1): 4 x0p - 4 x0n - 5 x1 - x2p + x2n >= -0
 constr(2): 3 x0p - 3 x0n - x1 + 3 x2p - 3 x2n + 2 x3 + 3 x4p - 3 x4n = 3
 constr(3): x0p - x0n + 5 x1 - 2 x2p + 2 x2n - x3 + 2 x4p - 2 x4n <= 21
 constr(4): - 5 x0p + 5 x0n + 5 x1 - 2 x2p + 2 x2n + 3 x3 >= -0
Bounds
 x0p <= 2
 x0n <= 3
 x1 <= 5
 x2p <= 1
 x2n <= 3
 x3 <= 5
 x4p <= 5
 x4n <= 2
Generals
 x0p x0n x1 x2p x2n x3 x4p x4n
End
