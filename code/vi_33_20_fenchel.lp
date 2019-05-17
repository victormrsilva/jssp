\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.4166666666666677 lambda(3,3,14) - 0.333333333333334 lambda(3,1,15)
   - 0.0833333333333313 lambda(2,3,17) - 0.166666666666667 lambda(1,1,14)
   - 0.583333333333335 lambda(2,3,15) - 0.0833333333333337 lambda(3,1,17)
Subject To
 sol(0): lambda(3,3,14) + lambda(3,1,15) + lambda(2,3,17) <= 1
 sol(1): lambda(3,3,14) + lambda(1,1,14) + lambda(2,3,15) <= 1
 sol(2): lambda(3,3,14) + lambda(2,3,17) + lambda(1,1,14) <= 1
 sol(3): lambda(3,3,14) + lambda(2,3,17) + lambda(3,1,17) <= 1
 sol(4): lambda(3,3,14) + lambda(3,1,15) + lambda(2,3,15) <= 1
 sol(5): lambda(3,3,14) + lambda(2,3,15) + lambda(3,1,17) <= 1
Bounds
 lambda(3,3,14) <= 1
 lambda(3,1,15) <= 1
 lambda(2,3,17) <= 1
 lambda(1,1,14) <= 1
 lambda(2,3,15) <= 1
 lambda(3,1,17) <= 1
End
