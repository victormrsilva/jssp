\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 5 x0p + 5 x0n - 5 x1p + 5 x1n + 5 x2p - 5 x2n + 3 x3p - 3 x3n
Subject To
 constr(0): x0p - x0n + 5 x2p - 5 x2n - 2 x3p + 2 x3n >= -11
 constr(1): - 2 x0p + 2 x0n - 4 x1p + 4 x1n + x2p - x2n - 2 x3p + 2 x3n
   >= -0
 constr(2): - 4 x0p + 4 x0n - 3 x1p + 3 x1n + 3 x2p - 3 x2n + 5 x3p - 5 x3n
   >= -13
 constr(3): - 4 x0p + 4 x0n + 4 x1p - 4 x1n - 2 x2p + 2 x2n - 3 x3p + 3 x3n
   <= 18
 constr(4): 4 x0p - 4 x0n - 3 x1p + 3 x1n + 3 x2p - 3 x2n + 4 x3p - 4 x3n
   <= 9
Bounds
 x0p <= 1
 x0n <= 3
 x1p <= 3
 x1n <= 4
 x2p <= 5
 x2n <= 3
 x3p <= 3
 x3n <= 5
End
