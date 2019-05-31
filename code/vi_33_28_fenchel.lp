\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.3333333333333335 lambda(1,1,27) - 0.3333333333333335 lambda(1,3,20)
   - 0.6666666666666667 lambda(3,1,20) - 0.3333333333333333 lambda(3,3,20)
   - 0.3333333333333333 lambda(3,1,26) - 0.3333333333333333 lambda(2,3,21)
Subject To
 sol(0): lambda(1,1,27) <= 1
 sol(1): lambda(1,1,27) + lambda(1,3,20) + lambda(3,1,20) <= 1
 sol(2): lambda(1,1,27) + lambda(3,3,20) <= 1
 sol(3): lambda(1,3,20) + lambda(3,1,26) <= 1
 sol(4): lambda(1,1,27) + lambda(3,1,20) <= 1
 sol(5): lambda(3,1,26) <= 1
 sol(6): lambda(1,1,27) + lambda(1,3,20) <= 1
 sol(7): lambda(1,1,27) + lambda(3,1,20) + lambda(2,3,21) <= 1
 sol(8): lambda(1,1,27) + lambda(2,3,21) <= 1
 sol(9): lambda(3,3,20) + lambda(3,1,26) <= 1
 sol(10): lambda(3,1,26) + lambda(2,3,21) <= 1
Bounds
 lambda(1,1,27) <= 1
 lambda(1,3,20) <= 1
 lambda(3,1,20) <= 1
 lambda(3,3,20) <= 1
 lambda(3,1,26) <= 1
 lambda(2,3,21) <= 1
End
