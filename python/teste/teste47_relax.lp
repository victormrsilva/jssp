\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 3 x0 + 4 x1p - 4 x1n - 5 x2p + 5 x2n + 4 x3
Subject To
 constr(0): 5 x0 + 4 x1p - 4 x1n + 4 x2p - 4 x2n - 4 x3 - 5 x4p + 5 x4n
   <= 17
 constr(1): x2p - x2n + 2 x3 - 5 x4p + 5 x4n <= 7
 constr(2): - 4 x0 + 5 x2p - 5 x2n - x3 + 5 x4p - 5 x4n = 8
Bounds
 x0 <= 2
 x1p <= 4
 x1n <= 1
 x2p <= 3
 x2n <= 4
 x3 <= 5
 x4p = 0
 x4n <= 2
End
