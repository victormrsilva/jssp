\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.125 lambda(3,1,24) - 0.375 lambda(2,3,24) - 0.375 lambda(1,1,24)
Subject To
 sol(0): lambda(3,1,24) <= 1
 sol(1): lambda(3,1,24) + lambda(2,3,24) <= 1
 sol(2): lambda(2,3,24) + lambda(1,1,24) <= 1
Bounds
 lambda(3,1,24) <= 1
 lambda(2,3,24) <= 1
 lambda(1,1,24) <= 1
End
