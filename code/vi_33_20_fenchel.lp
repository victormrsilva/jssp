\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.4798387096774193 lambda(1,2,0) - 0.1431451612903226 lambda(1,2,3)
   - 0.1028225806451613 lambda(3,3,3) - 0.1028225806451613 lambda(3,1,4)
   - 0.0403225806451613 lambda(1,2,6) - 0.1431451612903226 lambda(3,1,7)
Subject To
 sol(0): lambda(1,2,0) + lambda(1,2,3) + lambda(3,3,3) + lambda(3,1,4)
   + lambda(1,2,6) + lambda(3,1,7) <= 1
Bounds
 lambda(1,2,0) <= 1
 lambda(1,2,3) <= 1
 lambda(3,3,3) <= 1
 lambda(3,1,4) <= 1
 lambda(1,2,6) <= 1
 lambda(3,1,7) <= 1
End
