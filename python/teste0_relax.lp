\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0 + 5 x1 - x2p + x2n
Subject To
 constr(0): - 3 x1 - 5 x2p + 5 x2n >= -3
 constr(1): - 3 x0 - 5 x1 + 4 x2p - 4 x2n <= 3
 constr(2): 2 x0 + 2 x1 - 5 x2p + 5 x2n = 2
Bounds
 x0 <= 4
 x1 <= 1
 x2p <= 2
 x2n <= 2
End
