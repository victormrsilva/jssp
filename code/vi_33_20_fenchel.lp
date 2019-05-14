\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.1028225806451614 lambda(3,2,0) - 0.1431451612903226 lambda(3,2,3)
   - 0.4173387096774192 lambda(2,1,5) - 0.1431451612903226 lambda(1,3,6)
   - 0.0806451612903224 lambda(3,2,5) - 0.4798387096774194 lambda(1,3,3)
   - 0.4173387096774192 lambda(2,2,0) - 0.4798387096774194 lambda(1,2,0)
   - 0.1028225806451614 lambda(3,3,3) - 0.2963709677419356 lambda(2,2,3)
   - 0.1431451612903226 lambda(3,3,6) - 0.1431451612903226 lambda(3,1,7)
   - 0.1028225806451614 lambda(3,1,4) - 0.1431451612903226 lambda(1,2,3)
   - 0.2459677419354839 lambda(2,2,7) - 0.0403225806451612 lambda(1,2,6)
   - 0.3366935483870968 lambda(1,2,5)
Subject To
 sol(0): lambda(3,2,0) + lambda(3,2,3) + lambda(2,1,5) + lambda(1,3,6)
   <= 1
 sol(1): lambda(3,2,0) + lambda(2,1,5) + lambda(3,2,5) + lambda(1,3,3)
   <= 1
 sol(2): lambda(2,1,5) + lambda(3,2,5) + lambda(1,3,3) + lambda(2,2,0)
   <= 1
 sol(3): lambda(3,2,3) + lambda(2,1,5) + lambda(1,3,3) + lambda(1,2,0)
   <= 1
 sol(4): lambda(3,2,0) + lambda(2,1,5) + lambda(1,3,6) + lambda(3,2,5)
   + lambda(3,3,3) <= 1
 sol(5): lambda(1,2,0) + lambda(3,3,3) + lambda(2,2,3) + lambda(3,3,6)
   + lambda(3,1,7) <= 1
 sol(6): lambda(3,2,0) + lambda(1,3,6) + lambda(3,3,3) + lambda(2,2,3)
   + lambda(3,1,7) + lambda(3,1,4) <= 1
 sol(7): lambda(2,1,5) + lambda(1,3,6) + lambda(3,2,5) + lambda(2,2,0)
   + lambda(3,3,3) <= 1
 sol(8): lambda(3,2,0) + lambda(3,3,3) + lambda(2,2,3) + lambda(3,3,6)
   + lambda(3,1,7) <= 1
 sol(9): lambda(3,2,0) + lambda(1,3,6) + lambda(3,3,3) + lambda(3,1,7)
   + lambda(3,1,4) + lambda(1,2,3) + lambda(2,2,7) <= 1
 sol(10): lambda(3,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,3) + lambda(2,2,7) <= 1
 sol(11): lambda(3,2,0) + lambda(1,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(2,2,7) <= 1
 sol(12): lambda(3,2,3) + lambda(1,3,6) + lambda(1,2,0) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(13): lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,3) + lambda(3,1,7)
   + lambda(3,1,4) <= 1
 sol(14): lambda(3,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(15): lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,5) <= 1
 sol(16): lambda(3,2,0) + lambda(3,2,3) + lambda(2,1,5) + lambda(1,3,3)
   <= 1
 sol(17): lambda(1,3,6) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,1,7)
   + lambda(3,1,4) + lambda(1,2,3) + lambda(2,2,7) <= 1
 sol(18): lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(19): lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,3) + lambda(2,2,7) <= 1
 sol(20): lambda(3,2,3) + lambda(1,3,3) + lambda(1,2,0) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(21): lambda(3,2,3) + lambda(1,2,0) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,6) <= 1
 sol(22): lambda(2,1,5) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(23): lambda(1,3,3) + lambda(1,2,0) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(2,2,7) <= 1
 sol(24): lambda(3,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(25): lambda(1,3,6) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,2,3)
   + lambda(3,1,7) + lambda(3,1,4) <= 1
 sol(26): lambda(1,3,6) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,1,7)
   + lambda(3,1,4) + lambda(2,2,7) <= 1
 sol(27): lambda(2,1,5) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,5) <= 1
 sol(28): lambda(3,2,0) + lambda(2,1,5) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(29): lambda(1,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,3) + lambda(1,2,6) <= 1
 sol(30): lambda(3,2,0) + lambda(1,3,3) + lambda(2,2,3) + lambda(3,1,7)
   + lambda(3,1,4) <= 1
 sol(31): lambda(2,1,5) + lambda(1,3,6) + lambda(3,2,5) + lambda(1,2,0)
   + lambda(3,3,3) <= 1
 sol(32): lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(2,2,7) <= 1
 sol(33): lambda(2,1,5) + lambda(1,3,6) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(1,2,3) <= 1
 sol(34): lambda(2,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,5) <= 1
 sol(35): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(36): lambda(2,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,6) <= 1
 sol(37): lambda(3,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,5) <= 1
 sol(38): lambda(3,2,0) + lambda(3,2,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,6) <= 1
 sol(39): lambda(3,2,0) + lambda(2,1,5) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,5) <= 1
 sol(40): lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,5) <= 1
 sol(41): lambda(2,1,5) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,6) <= 1
 sol(42): lambda(2,1,5) + lambda(3,2,5) + lambda(1,3,3) + lambda(1,2,0)
   <= 1
 sol(43): lambda(2,1,5) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,5) <= 1
 sol(44): lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(1,2,6) <= 1
 sol(45): lambda(3,2,0) + lambda(3,2,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(1,2,6) <= 1
 sol(46): lambda(1,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,5) <= 1
 sol(47): lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(48): lambda(3,2,0) + lambda(2,1,5) + lambda(1,3,6) + lambda(3,3,3)
   + lambda(1,2,3) <= 1
 sol(49): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(50): lambda(3,2,3) + lambda(1,2,0) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(51): lambda(3,2,0) + lambda(3,2,3) + lambda(3,3,6) + lambda(3,1,7)
   + lambda(2,2,7) <= 1
 sol(52): lambda(3,2,3) + lambda(2,1,5) + lambda(1,3,6) + lambda(1,2,0)
   <= 1
 sol(53): lambda(3,2,3) + lambda(2,1,5) + lambda(1,2,0) + lambda(3,3,6)
   + lambda(1,2,6) <= 1
 sol(54): lambda(3,2,0) + lambda(3,3,3) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,5) <= 1
Bounds
 lambda(3,2,0) <= 1
 lambda(3,2,3) <= 1
 lambda(2,1,5) <= 1
 lambda(1,3,6) <= 1
 lambda(3,2,5) <= 1
 lambda(1,3,3) <= 1
 lambda(2,2,0) <= 1
 lambda(1,2,0) <= 1
 lambda(3,3,3) <= 1
 lambda(2,2,3) <= 1
 lambda(3,3,6) <= 1
 lambda(3,1,7) <= 1
 lambda(3,1,4) <= 1
 lambda(1,2,3) <= 1
 lambda(2,2,7) <= 1
 lambda(1,2,6) <= 1
 lambda(1,2,5) <= 1
End
