\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 4 x0p + 4 x0n + x1p - x1n + 2 x2 + 2 x3p - 2 x3n
Subject To
 constr(0): - 4 x0p + 4 x0n - 3 x1p + 3 x1n + 4 x2 - 5 x3p + 5 x3n - 3 x4p
   + 3 x4n >= -1
 constr(1): 4 x0p - 4 x0n - 4 x1p + 4 x1n - 5 x2 - 3 x4p + 3 x4n <= -0
 constr(2): 2 x0p - 2 x0n - 2 x1p + 2 x1n - 5 x3p + 5 x3n - 3 x4p + 3 x4n
   <= 3
 constr(3): 3 x0p - 3 x0n + 4 x1p - 4 x1n + 5 x2 + 3 x3p - 3 x3n + 2 x4p
   - 2 x4n >= -0
 constr(4): - x1p + x1n + 4 x2 - 2 x3p + 2 x3n = 2
Bounds
 x0p = 0
 x0n <= 4
 x1p = 0
 x1n <= 1
 x2 <= 1
 x3p <= 1
 x3n <= 3
 x4p <= 1
 x4n <= 4
End
