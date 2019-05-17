\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.2777777777777778 lambda(1,3,11) - 0.0555555555555553 lambda(3,2,13)
   - 0.1666666666666665 lambda(3,3,10) - 0.3333333333333331 lambda(2,2,12)
   - 0.1111111111111112 lambda(1,1,14) - 0.388888888888889 lambda(3,2,9)
   - 0.4444444444444443 lambda(2,3,13) - 0.2777777777777777 lambda(3,1,14)
   - 0.5555555555555557 lambda(1,1,9) - 0.1111111111111113 lambda(3,2,12)
   - 0.1111111111111113 lambda(3,2,17) - 0.0555555555555553 lambda(1,2,9)
   - 0.0555555555555552 lambda(3,2,14) - 0.0555555555555554 lambda(1,1,13)
   - 0.2222222222222226 lambda(3,1,17) - 0.2222222222222225 lambda(2,2,14)
   - 0.2222222222222225 lambda(3,2,15) - 0.1111111111111112 lambda(1,2,10)
Subject To
 sol(0): lambda(1,3,11) + lambda(3,2,13) <= 1
 sol(1): lambda(3,3,10) + lambda(2,2,12) + lambda(1,1,14) <= 1
 sol(2): lambda(3,2,9) + lambda(2,3,13) + lambda(3,1,14) <= 1
 sol(3): lambda(3,2,13) + lambda(2,3,13) + lambda(1,1,9) <= 1
 sol(4): lambda(1,3,11) + lambda(3,2,12) + lambda(3,2,17) <= 1
 sol(5): lambda(2,3,13) + lambda(1,2,9) + lambda(3,2,14) + lambda(1,1,13)
   <= 1
 sol(6): lambda(2,3,13) + lambda(1,1,9) + lambda(3,2,12) + lambda(3,1,17)
   <= 1
 sol(7): lambda(1,1,14) + lambda(2,3,13) + lambda(1,2,9) + lambda(3,2,14)
   <= 1
 sol(8): lambda(3,2,13) + lambda(2,3,13) + lambda(1,2,9) + lambda(1,1,13)
   <= 1
 sol(9): lambda(1,3,11) + lambda(3,2,9) + lambda(3,2,14) <= 1
 sol(10): lambda(3,3,10) + lambda(1,2,9) + lambda(3,1,17) + lambda(2,2,14)
   <= 1
 sol(11): lambda(3,2,9) + lambda(2,3,13) + lambda(1,1,9) + lambda(3,2,17)
   <= 1
 sol(12): lambda(3,3,10) + lambda(1,1,14) + lambda(3,2,17) + lambda(1,2,9)
   <= 1
 sol(13): lambda(3,3,10) + lambda(1,1,9) + lambda(3,1,17) + lambda(2,2,14)
   <= 1
 sol(14): lambda(3,3,10) + lambda(2,2,12) + lambda(1,1,9) + lambda(3,1,17)
   <= 1
 sol(15): lambda(1,3,11) + lambda(3,2,12) + lambda(3,1,17) <= 1
 sol(16): lambda(3,3,10) + lambda(1,1,14) + lambda(1,2,9) + lambda(2,2,14)
   <= 1
 sol(17): lambda(1,1,14) + lambda(2,3,13) + lambda(1,2,9) + lambda(3,2,15)
   <= 1
 sol(18): lambda(3,3,10) + lambda(3,2,17) + lambda(1,2,9) + lambda(1,1,13)
   <= 1
 sol(19): lambda(1,1,14) + lambda(2,3,13) + lambda(3,2,17) + lambda(1,2,9)
   <= 1
 sol(20): lambda(3,2,13) + lambda(1,1,14) + lambda(2,3,13) + lambda(1,2,9)
   <= 1
 sol(21): lambda(3,2,9) + lambda(2,3,13) + lambda(1,1,9) + lambda(3,2,15)
   <= 1
 sol(22): lambda(1,3,11) + lambda(3,2,9) + lambda(3,1,14) + lambda(2,2,14)
   <= 1
 sol(23): lambda(1,3,11) + lambda(3,2,9) + lambda(3,2,15) <= 1
 sol(24): lambda(3,2,9) + lambda(2,3,13) + lambda(1,1,9) + lambda(3,2,14)
   <= 1
 sol(25): lambda(2,3,13) + lambda(3,1,14) + lambda(1,2,9) <= 1
 sol(26): lambda(3,2,9) + lambda(2,3,13) + lambda(1,1,9) + lambda(3,1,17)
   <= 1
 sol(27): lambda(1,3,11) + lambda(2,2,12) + lambda(3,1,14) <= 1
 sol(28): lambda(2,3,13) + lambda(1,1,9) + lambda(3,2,12) + lambda(3,2,17)
   <= 1
 sol(29): lambda(2,3,13) + lambda(3,1,14) + lambda(1,2,10) <= 1
 sol(30): lambda(3,3,10) + lambda(1,2,9) + lambda(1,1,13) + lambda(2,2,14)
   <= 1
 sol(31): lambda(3,3,10) + lambda(1,1,14) + lambda(3,2,17) + lambda(1,2,10)
   <= 1
 sol(32): lambda(3,1,14) + lambda(2,2,14) + lambda(1,2,10) <= 1
 sol(33): lambda(3,2,9) + lambda(1,1,13) + lambda(2,2,14) <= 1
 sol(34): lambda(1,3,11) + lambda(3,2,9) + lambda(3,1,17) + lambda(2,2,14)
   <= 1
 sol(35): lambda(2,3,13) + lambda(1,2,9) + lambda(1,1,13) + lambda(3,2,15)
   <= 1
 sol(36): lambda(3,2,9) + lambda(2,3,13) + lambda(1,1,13) + lambda(3,2,15)
   <= 1
 sol(37): lambda(3,3,10) + lambda(2,2,12) + lambda(1,1,13) <= 1
 sol(38): lambda(1,1,14) + lambda(3,2,9) + lambda(2,2,14) <= 1
 sol(39): lambda(1,3,11) + lambda(3,2,9) + lambda(3,2,17) <= 1
 sol(40): lambda(3,2,9) + lambda(2,3,13) + lambda(3,2,17) + lambda(1,1,13)
   <= 1
 sol(41): lambda(2,3,13) + lambda(3,1,17) + lambda(1,2,10) <= 1
 sol(42): lambda(1,1,14) + lambda(2,3,13) + lambda(3,2,14) + lambda(1,2,10)
   <= 1
 sol(43): lambda(3,3,10) + lambda(1,1,9) + lambda(3,2,17) <= 1
 sol(44): lambda(3,2,9) + lambda(2,3,13) + lambda(3,2,14) + lambda(1,1,13)
   <= 1
 sol(45): lambda(1,1,14) + lambda(3,2,9) + lambda(2,3,13) + lambda(3,2,14)
   <= 1
 sol(46): lambda(3,1,14) + lambda(1,2,9) + lambda(2,2,14) <= 1
 sol(47): lambda(3,2,9) + lambda(1,1,9) + lambda(3,1,17) + lambda(2,2,14)
   <= 1
 sol(48): lambda(1,1,14) + lambda(3,2,9) + lambda(2,3,13) + lambda(3,2,17)
   <= 1
 sol(49): lambda(3,3,10) + lambda(1,1,14) + lambda(2,2,14) + lambda(1,2,10)
   <= 1
 sol(50): lambda(2,3,13) + lambda(3,2,17) + lambda(1,2,9) + lambda(1,1,13)
   <= 1
 sol(51): lambda(1,1,14) + lambda(2,3,13) + lambda(3,2,15) + lambda(1,2,10)
   <= 1
 sol(52): lambda(3,3,10) + lambda(3,1,17) + lambda(2,2,14) + lambda(1,2,10)
   <= 1
 sol(53): lambda(1,3,11) + lambda(2,2,12) + lambda(3,1,17) <= 1
 sol(54): lambda(1,1,14) + lambda(2,3,13) + lambda(3,2,17) + lambda(1,2,10)
   <= 1
 sol(55): lambda(1,1,14) + lambda(3,2,9) + lambda(2,3,13) + lambda(3,2,15)
   <= 1
 sol(56): lambda(1,1,14) + lambda(2,3,13) + lambda(3,2,12) + lambda(3,2,17)
   <= 1
 sol(57): lambda(2,3,13) + lambda(3,2,12) + lambda(3,2,17) + lambda(1,1,13)
   <= 1
 sol(58): lambda(2,3,13) + lambda(1,2,9) + lambda(3,1,17) <= 1
Bounds
 lambda(1,3,11) <= 1
 lambda(3,2,13) <= 1
 lambda(3,3,10) <= 1
 lambda(2,2,12) <= 1
 lambda(1,1,14) <= 1
 lambda(3,2,9) <= 1
 lambda(2,3,13) <= 1
 lambda(3,1,14) <= 1
 lambda(1,1,9) <= 1
 lambda(3,2,12) <= 1
 lambda(3,2,17) <= 1
 lambda(1,2,9) <= 1
 lambda(3,2,14) <= 1
 lambda(1,1,13) <= 1
 lambda(3,1,17) <= 1
 lambda(2,2,14) <= 1
 lambda(3,2,15) <= 1
 lambda(1,2,10) <= 1
End
