\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  5 x0p - 5 x0n - 2 x3p + 2 x3n
Subject To
 constr(0): 4 x0p - 4 x0n + 4 x1 - 3 x2p + 3 x2n - 2 x3p + 2 x3n >= -19
 constr(1): - 5 x0p + 5 x0n - 3 x2p + 3 x2n - 3 x3p + 3 x3n = 18
 constr(2): - 5 x0p + 5 x0n - 5 x1 + 3 x2p - 3 x2n <= 8
 constr(3): - 4 x0p + 4 x0n + 3 x1 - 3 x2p + 3 x2n - 3 x3p + 3 x3n <= 15
 constr(4): 5 x0p - 5 x0n + 3 x1 + 3 x2p - 3 x2n + 2 x3p - 2 x3n <= 11
Bounds
 x0p <= 1
 x0n <= 4
 x1 <= 1
 x2p <= 3
 x2n <= 4
 x3p <= 2
 x3n <= 2
Generals
 x0p x0n x1 x2p x2n x3p x3n
End
