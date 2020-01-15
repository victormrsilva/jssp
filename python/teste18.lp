\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - 5 x1 + x2p - x2n - x3p + x3n
Subject To
 constr(0): - x1 + 5 x2p - 5 x2n + 3 x3p - 3 x3n >= -12
 constr(1): - 4 x0p + 4 x0n + 3 x1 - x2p + x2n - 2 x3p + 2 x3n <= 6
 constr(2): - 3 x0p + 3 x0n + 4 x1 + 2 x2p - 2 x2n - 4 x3p + 4 x3n = 10
 constr(3): 5 x0p - 5 x0n + 4 x1 + 3 x2p - 3 x2n - x3p + x3n >= -9
 constr(4): 4 x0p - 4 x0n - 2 x1 + 5 x2p - 5 x2n - 5 x3p + 5 x3n <= 12
Bounds
 x0p <= 2
 x0n <= 1
 x1 <= 1
 x2p <= 4
 x2n <= 3
 x3p <= 4
 x3n <= 5
Generals
 x0p x0n x1 x2p x2n x3p x3n
End
