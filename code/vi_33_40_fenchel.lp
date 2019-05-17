\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.1666666666666685 lambda(2,1,27) - 0.3333333333333315 lambda(3,3,29)
   - 0.1666666666666685 lambda(3,2,21) - 0.166666666666663 lambda(1,1,22)
   - 0.3333333333333315 lambda(2,3,26) - 0.1666666666666685 lambda(3,3,30)
   - 0.1666666666666685 lambda(2,3,34) - 0.166666666666663 lambda(1,3,30)
   - 0.1666666666666685 lambda(1,3,24) - 0.1666666666666685 lambda(3,3,26)
Subject To
 sol(0): lambda(2,1,27) + lambda(3,3,29) <= 1
 sol(1): lambda(3,2,21) + lambda(1,1,22) + lambda(2,3,26) + lambda(3,3,30)
   <= 1
 sol(2): lambda(3,2,21) + lambda(1,1,22) + lambda(2,3,26) + lambda(2,3,34)
   <= 1
 sol(3): lambda(3,2,21) + lambda(1,1,22) + lambda(2,3,26) + lambda(1,3,30)
   <= 1
 sol(4): lambda(3,2,21) + lambda(2,3,34) + lambda(1,3,24) <= 1
 sol(5): lambda(2,1,27) + lambda(3,2,21) + lambda(1,3,30) <= 1
 sol(6): lambda(1,1,22) + lambda(3,3,26) <= 1
 sol(7): lambda(3,3,29) + lambda(1,1,22) <= 1
 sol(8): lambda(2,1,27) + lambda(3,3,26) <= 1
 sol(9): lambda(2,1,27) + lambda(3,2,21) + lambda(1,3,24) <= 1
 sol(10): lambda(2,1,27) + lambda(3,2,21) + lambda(3,3,30) <= 1
Bounds
 lambda(2,1,27) <= 1
 lambda(3,3,29) <= 1
 lambda(3,2,21) <= 1
 lambda(1,1,22) <= 1
 lambda(2,3,26) <= 1
 lambda(3,3,30) <= 1
 lambda(2,3,34) <= 1
 lambda(1,3,30) <= 1
 lambda(1,3,24) <= 1
 lambda(3,3,26) <= 1
End
