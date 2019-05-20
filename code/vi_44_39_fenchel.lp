\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.333333333333333 lambda(3,1,37) - 0.3333333333333356 lambda(3,1,31)
   - 0.3333333333333329 lambda(2,3,31) - 0.6666666666666684 lambda(1,2,31)
   - 0.333333333333333 lambda(3,2,30) - 0.6666666666666663 lambda(4,3,29)
   - 0.3333333333333329 lambda(2,2,26) - 0.3333333333333316 lambda(3,3,26)
   - 0.3333333333333328 lambda(3,1,26)
Subject To
 sol(0): lambda(3,1,37) <= 1
 sol(1): lambda(3,1,31) <= 1
 sol(2): lambda(3,1,31) + lambda(2,3,31) <= 1
 sol(3): lambda(3,1,31) + lambda(1,2,31) <= 1
 sol(4): lambda(3,1,37) + lambda(2,3,31) + lambda(1,2,31) <= 1
 sol(5): lambda(3,1,37) + lambda(3,2,30) <= 1
 sol(6): lambda(3,1,37) + lambda(2,3,31) + lambda(3,2,30) <= 1
 sol(7): lambda(3,1,31) + lambda(4,3,29) <= 1
 sol(8): lambda(3,1,37) + lambda(1,2,31) + lambda(4,3,29) <= 1
 sol(9): lambda(3,1,31) + lambda(1,2,31) + lambda(4,3,29) <= 1
 sol(10): lambda(3,1,37) + lambda(2,3,31) + lambda(1,2,31) + lambda(2,2,26)
   <= 1
 sol(11): lambda(3,1,31) + lambda(4,3,29) + lambda(2,2,26) <= 1
 sol(12): lambda(3,1,31) + lambda(2,3,31) + lambda(1,2,31) <= 1
 sol(13): lambda(3,1,37) + lambda(1,2,31) + lambda(4,3,29) + lambda(2,2,26)
   <= 1
 sol(14): lambda(3,1,31) + lambda(2,2,26) + lambda(3,3,26) <= 1
 sol(15): lambda(3,1,37) + lambda(1,2,31) + lambda(4,3,29) + lambda(2,2,26)
   + lambda(3,3,26) <= 1
 sol(16): lambda(3,1,37) + lambda(2,3,31) <= 1
 sol(17): lambda(2,3,31) + lambda(2,2,26) + lambda(3,1,26) <= 1
 sol(18): lambda(3,1,37) + lambda(1,2,31) <= 1
 sol(19): lambda(3,1,31) + lambda(1,2,31) + lambda(4,3,29) + lambda(2,2,26)
   + lambda(3,3,26) <= 1
 sol(20): lambda(3,1,37) + lambda(2,3,31) + lambda(1,2,31) + lambda(2,2,26)
   + lambda(3,3,26) <= 1
 sol(21): lambda(2,3,31) + lambda(1,2,31) + lambda(2,2,26) + lambda(3,1,26)
   <= 1
 sol(22): lambda(3,1,31) + lambda(1,2,31) + lambda(2,2,26) <= 1
 sol(23): lambda(3,1,31) + lambda(1,2,31) + lambda(4,3,29) + lambda(2,2,26)
   <= 1
 sol(24): lambda(3,1,31) + lambda(2,3,31) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(25): lambda(3,1,37) + lambda(3,2,30) + lambda(4,3,29) <= 1
 sol(26): lambda(3,1,37) + lambda(2,2,26) + lambda(3,3,26) <= 1
 sol(27): lambda(3,1,37) + lambda(4,3,29) + lambda(2,2,26) <= 1
 sol(28): lambda(3,1,37) + lambda(3,3,26) <= 1
 sol(29): lambda(3,1,31) + lambda(1,2,31) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(30): lambda(3,1,37) + lambda(4,3,29) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(31): lambda(3,1,37) + lambda(2,3,31) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(32): lambda(3,1,37) + lambda(1,2,31) + lambda(2,2,26) <= 1
 sol(33): lambda(3,1,31) + lambda(4,3,29) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(34): lambda(3,1,31) + lambda(2,3,31) + lambda(2,2,26) <= 1
 sol(35): lambda(3,1,37) + lambda(1,2,31) + lambda(2,2,26) + lambda(3,3,26)
   <= 1
 sol(36): lambda(3,1,31) + lambda(2,2,26) <= 1
 sol(37): lambda(3,1,37) + lambda(2,3,31) + lambda(2,2,26) <= 1
 sol(38): lambda(3,1,37) + lambda(2,2,26) <= 1
 sol(39): lambda(1,2,31) + lambda(4,3,29) + lambda(3,1,26) <= 1
 sol(40): lambda(3,1,37) + lambda(4,3,29) <= 1
 sol(41): lambda(3,1,31) + lambda(2,3,31) + lambda(1,2,31) + lambda(2,2,26)
   <= 1
 sol(42): lambda(3,1,31) + lambda(4,3,29) + lambda(3,3,26) <= 1
 sol(43): lambda(1,2,31) + lambda(4,3,29) + lambda(2,2,26) + lambda(3,1,26)
   <= 1
 sol(44): lambda(2,3,31) + lambda(1,2,31) + lambda(3,1,26) <= 1
 sol(45): lambda(3,1,37) + lambda(3,2,30) + lambda(3,3,26) <= 1
 sol(46): lambda(3,1,31) + lambda(2,3,31) + lambda(1,2,31) + lambda(2,2,26)
   + lambda(3,3,26) <= 1
 sol(47): lambda(3,1,31) + lambda(1,2,31) + lambda(4,3,29) + lambda(3,3,26)
   <= 1
 sol(48): lambda(3,1,31) + lambda(2,3,31) + lambda(3,3,26) <= 1
 sol(49): lambda(3,1,37) + lambda(3,2,30) + lambda(4,3,29) + lambda(3,3,26)
   <= 1
 sol(50): lambda(3,1,37) + lambda(1,2,31) + lambda(4,3,29) + lambda(3,3,26)
   <= 1
 sol(51): lambda(3,1,37) + lambda(4,3,29) + lambda(3,3,26) <= 1
 sol(52): lambda(2,3,31) + lambda(3,1,26) <= 1
 sol(53): lambda(3,1,31) + lambda(3,3,26) <= 1
 sol(54): lambda(3,1,37) + lambda(2,3,31) + lambda(3,3,26) <= 1
 sol(55): lambda(3,1,37) + lambda(2,3,31) + lambda(3,2,30) + lambda(3,3,26)
   <= 1
 sol(56): lambda(3,1,31) + lambda(2,3,31) + lambda(1,2,31) + lambda(3,3,26)
   <= 1
 sol(57): lambda(3,1,37) + lambda(2,3,31) + lambda(1,2,31) + lambda(3,3,26)
   <= 1
 sol(58): lambda(3,1,31) + lambda(1,2,31) + lambda(3,3,26) <= 1
 sol(59): lambda(3,1,37) + lambda(1,2,31) + lambda(3,3,26) <= 1
Bounds
 lambda(3,1,37) <= 1
 lambda(3,1,31) <= 1
 lambda(2,3,31) <= 1
 lambda(1,2,31) <= 1
 lambda(3,2,30) <= 1
 lambda(4,3,29) <= 1
 lambda(2,2,26) <= 1
 lambda(3,3,26) <= 1
 lambda(3,1,26) <= 1
End
