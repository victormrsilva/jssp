\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.1666666666666667 lambda(3,1,17) - 0.1666666666666673 lambda(3,1,16)
   - 0.1666666666666667 lambda(2,3,17) - 0.333333333333334 lambda(3,3,15)
   - 0.6666666666666666 lambda(2,3,16) - 0.166666666666666 lambda(3,3,14)
   - 0.166666666666666 lambda(3,1,15) - 0.1666666666666664 lambda(1,1,15)
   - 0.1666666666666667 lambda(1,1,14)
Subject To
 sol(0): lambda(3,1,17) <= 1
 sol(1): lambda(3,1,16) + lambda(2,3,17) <= 1
 sol(2): lambda(3,1,17) + lambda(3,3,15) <= 1
 sol(3): lambda(3,1,16) + lambda(2,3,16) <= 1
 sol(4): lambda(3,1,17) + lambda(2,3,17) + lambda(3,3,15) <= 1
 sol(5): lambda(3,1,17) + lambda(3,3,15) + lambda(2,3,16) <= 1
 sol(6): lambda(2,3,17) + lambda(3,3,14) + lambda(3,1,15) <= 1
 sol(7): lambda(2,3,16) + lambda(1,1,15) <= 1
 sol(8): lambda(2,3,16) + lambda(3,3,14) + lambda(3,1,15) <= 1
 sol(9): lambda(2,3,17) + lambda(1,1,15) <= 1
 sol(10): lambda(3,1,17) + lambda(2,3,17) <= 1
 sol(11): lambda(3,1,16) + lambda(2,3,17) + lambda(3,3,15) <= 1
 sol(12): lambda(2,3,16) + lambda(3,1,15) <= 1
 sol(13): lambda(2,3,17) + lambda(3,3,14) + lambda(1,1,15) <= 1
 sol(14): lambda(2,3,17) + lambda(1,1,14) <= 1
 sol(15): lambda(2,3,16) + lambda(1,1,14) <= 1
 sol(16): lambda(3,1,17) + lambda(2,3,16) <= 1
 sol(17): lambda(2,3,16) + lambda(3,3,14) + lambda(1,1,15) <= 1
 sol(18): lambda(3,1,17) + lambda(2,3,16) + lambda(3,3,14) <= 1
 sol(19): lambda(2,3,17) + lambda(3,3,15) + lambda(1,1,15) <= 1
 sol(20): lambda(2,3,16) + lambda(3,3,14) + lambda(1,1,14) <= 1
 sol(21): lambda(2,3,17) + lambda(3,1,15) <= 1
 sol(22): lambda(3,3,15) + lambda(2,3,16) + lambda(1,1,14) <= 1
 sol(23): lambda(3,1,17) + lambda(2,3,17) + lambda(3,3,14) <= 1
 sol(24): lambda(2,3,17) + lambda(3,3,15) + lambda(1,1,14) <= 1
 sol(25): lambda(3,1,16) + lambda(2,3,16) + lambda(3,3,14) <= 1
 sol(26): lambda(3,1,16) + lambda(2,3,17) + lambda(3,3,14) <= 1
 sol(27): lambda(3,1,17) + lambda(3,3,14) <= 1
 sol(28): lambda(2,3,17) + lambda(3,3,14) + lambda(1,1,14) <= 1
 sol(29): lambda(3,1,16) + lambda(3,3,15) + lambda(2,3,16) <= 1
 sol(30): lambda(3,3,15) + lambda(2,3,16) + lambda(1,1,15) <= 1
Bounds
 lambda(3,1,17) <= 1
 lambda(3,1,16) <= 1
 lambda(2,3,17) <= 1
 lambda(3,3,15) <= 1
 lambda(2,3,16) <= 1
 lambda(3,3,14) <= 1
 lambda(3,1,15) <= 1
 lambda(1,1,15) <= 1
 lambda(1,1,14) <= 1
End
