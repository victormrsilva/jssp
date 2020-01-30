\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0p + 3 x0n + 4 x1 + 4 x3p - 4 x3n - 2 x4p + 2 x4n
Subject To
 constr(0): 5 x0p - 5 x0n + 3 x1 + 5 x2 - 2 x3p + 2 x3n <= 12
 constr(1): 5 x0p - 5 x0n + 2 x1 - 4 x2 + 5 x3p - 5 x3n - 5 x4p + 5 x4n
   >= -2
 constr(2): - 2 x0p + 2 x0n - 5 x1 + 3 x2 + 3 x3p - 3 x3n + x4p - x4n
   <= 15
Bounds
 x0p = 0
 x0n <= 4
 x1 <= 4
 x2 <= 4
 x3p <= 2
 x3n <= 3
 x4p <= 1
 x4n <= 4
Generals
 x0p x0n x1 x2 x3p x3n x4p x4n
End
