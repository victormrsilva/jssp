\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.4 lambda(2,2,19) - 0.4 lambda(1,3,19) - 0.4 lambda(2,3,24)
   - 0.2 lambda(2,1,18) - 0.2 lambda(3,1,21) - 0.4 lambda(1,1,21)
   - 0.2 lambda(3,1,19) - 0.4 lambda(2,1,25)
Subject To
 sol(0): lambda(2,2,19) + lambda(1,3,19) + lambda(2,3,24) + lambda(2,1,18)
   + lambda(3,1,21) <= 1
 sol(1): lambda(2,2,19) + lambda(1,3,19) + lambda(2,3,24) + lambda(2,1,18)
   + lambda(1,1,21) <= 1
 sol(2): lambda(2,2,19) + lambda(1,3,19) + lambda(2,3,24) + lambda(2,1,18)
   + lambda(3,1,19) + lambda(2,1,25) <= 1
Bounds
 lambda(2,2,19) <= 1
 lambda(1,3,19) <= 1
 lambda(2,3,24) <= 1
 lambda(2,1,18) <= 1
 lambda(3,1,21) <= 1
 lambda(1,1,21) <= 1
 lambda(3,1,19) <= 1
 lambda(2,1,25) <= 1
End
