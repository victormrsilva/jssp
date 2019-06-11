\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.3333333333333328 lambda(2,3,17) - 0.3333333333333315 lambda(3,1,17)
   - 0.3333333333333315 lambda(3,3,16) - 0.3333333333333324 lambda(1,1,15)
   - 0.3333333333333347 lambda(2,3,14) - 0.3333333333333315 lambda(3,2,13)
   - 0.1666666666666678 lambda(2,3,15) - 0.3333333333333328 lambda(2,1,13)
   - 0.1666666666666678 lambda(2,1,11) - 0.1666666666666655 lambda(3,3,11)
   - 0.1666666666666655 lambda(3,1,12) - 0.0416666666666636 lambda(3,1,10)
   - 0.1666666666666646 lambda(2,3,12) - 0.3333333333333348 lambda(1,1,13)
Subject To
 sol(0): lambda(2,3,17) <= 1
 sol(1): lambda(2,3,17) + lambda(3,1,17) <= 1
 sol(2): lambda(2,3,17) + lambda(3,3,16) <= 1
 sol(3): lambda(2,3,17) + lambda(1,1,15) <= 1
 sol(4): lambda(2,3,17) + lambda(3,3,16) + lambda(1,1,15) <= 1
 sol(5): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) <= 1
 sol(6): lambda(1,1,15) + lambda(2,3,14) <= 1
 sol(7): lambda(2,3,17) + lambda(3,2,13) <= 1
 sol(8): lambda(2,3,17) + lambda(3,3,16) + lambda(3,2,13) <= 1
 sol(9): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) + lambda(3,2,13)
   <= 1
 sol(10): lambda(3,1,17) + lambda(3,2,13) + lambda(2,3,15) <= 1
 sol(11): lambda(3,1,17) + lambda(2,3,14) + lambda(3,2,13) <= 1
 sol(12): lambda(1,1,15) + lambda(2,3,14) + lambda(3,2,13) <= 1
 sol(13): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) + lambda(2,1,13)
   <= 1
 sol(14): lambda(3,1,17) + lambda(2,3,14) <= 1
 sol(15): lambda(2,3,17) + lambda(3,2,13) + lambda(2,1,13) <= 1
 sol(16): lambda(2,3,17) + lambda(2,1,11) <= 1
 sol(17): lambda(1,1,15) + lambda(2,3,15) + lambda(2,1,11) <= 1
 sol(18): lambda(3,1,17) + lambda(2,3,15) + lambda(2,1,11) <= 1
 sol(19): lambda(1,1,15) + lambda(2,3,15) + lambda(3,3,11) + lambda(3,1,12)
   <= 1
 sol(20): lambda(3,1,10) + lambda(2,3,12) + lambda(1,1,13) <= 1
 sol(21): lambda(2,3,17) + lambda(3,3,16) + lambda(3,2,13) + lambda(1,1,13)
   <= 1
 sol(22): lambda(1,1,15) + lambda(2,3,14) + lambda(3,3,11) + lambda(3,1,12)
   <= 1
 sol(23): lambda(2,3,17) + lambda(1,1,15) + lambda(2,1,11) <= 1
 sol(24): lambda(2,3,17) + lambda(3,1,17) + lambda(2,1,11) <= 1
 sol(25): lambda(2,3,15) + lambda(1,1,13) <= 1
 sol(26): lambda(2,3,17) + lambda(3,1,17) + lambda(3,2,13) + lambda(2,1,11)
   <= 1
 sol(27): lambda(2,3,17) + lambda(3,3,16) + lambda(3,2,13) + lambda(2,1,11)
   <= 1
 sol(28): lambda(1,1,15) + lambda(3,2,13) + lambda(2,3,15) + lambda(2,1,11)
   <= 1
 sol(29): lambda(3,1,17) + lambda(2,3,15) <= 1
 sol(30): lambda(2,3,17) + lambda(2,1,13) + lambda(3,1,10) <= 1
 sol(31): lambda(1,1,15) + lambda(3,1,12) + lambda(2,3,12) <= 1
 sol(32): lambda(2,3,17) + lambda(1,1,13) <= 1
 sol(33): lambda(2,3,17) + lambda(1,1,15) + lambda(3,2,13) + lambda(2,1,11)
   <= 1
 sol(34): lambda(1,1,15) + lambda(2,3,15) + lambda(2,1,11) + lambda(3,3,11)
   <= 1
 sol(35): lambda(2,3,17) + lambda(3,3,16) + lambda(2,1,11) <= 1
 sol(36): lambda(3,3,11) + lambda(2,3,12) + lambda(1,1,13) <= 1
 sol(37): lambda(3,1,17) + lambda(2,3,12) <= 1
 sol(38): lambda(1,1,15) + lambda(3,3,11) + lambda(3,1,12) + lambda(2,3,12)
   <= 1
 sol(39): lambda(2,3,14) + lambda(3,1,10) + lambda(1,1,13) <= 1
 sol(40): lambda(2,3,17) + lambda(3,3,16) + lambda(1,1,15) + lambda(2,1,11)
   <= 1
 sol(41): lambda(2,3,17) + lambda(3,3,16) + lambda(1,1,15) + lambda(3,2,13)
   <= 1
 sol(42): lambda(3,1,17) + lambda(3,2,13) + lambda(2,3,15) + lambda(2,1,11)
   <= 1
 sol(43): lambda(1,1,15) + lambda(2,3,15) <= 1
 sol(44): lambda(2,3,17) + lambda(3,3,16) + lambda(2,1,13) <= 1
 sol(45): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) + lambda(3,2,13)
   + lambda(2,1,13) <= 1
 sol(46): lambda(2,3,17) + lambda(3,3,16) + lambda(3,2,13) + lambda(2,1,13)
   <= 1
 sol(47): lambda(2,3,17) + lambda(1,1,15) + lambda(3,3,11) + lambda(3,1,12)
   <= 1
 sol(48): lambda(2,3,15) + lambda(3,1,10) + lambda(1,1,13) <= 1
 sol(49): lambda(1,1,15) + lambda(2,3,15) + lambda(3,1,10) <= 1
 sol(50): lambda(2,3,17) + lambda(3,1,17) + lambda(2,1,13) + lambda(3,3,11)
   <= 1
 sol(51): lambda(1,1,15) + lambda(3,1,10) + lambda(2,3,12) <= 1
 sol(52): lambda(1,1,15) + lambda(2,3,14) + lambda(3,1,10) <= 1
 sol(53): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) + lambda(3,2,13)
   + lambda(2,1,11) <= 1
 sol(54): lambda(2,3,15) + lambda(3,1,10) <= 1
 sol(55): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,16) + lambda(2,1,11)
   <= 1
 sol(56): lambda(2,3,17) + lambda(1,1,15) + lambda(3,1,10) <= 1
 sol(57): lambda(2,3,17) + lambda(3,1,17) + lambda(2,1,13) <= 1
 sol(58): lambda(2,3,17) + lambda(3,1,12) <= 1
 sol(59): lambda(3,1,17) + lambda(2,3,15) + lambda(2,1,11) + lambda(3,3,11)
   <= 1
 sol(60): lambda(2,3,17) + lambda(3,1,17) + lambda(2,1,11) + lambda(3,3,11)
   <= 1
 sol(61): lambda(2,3,17) + lambda(1,1,15) + lambda(3,3,11) <= 1
 sol(62): lambda(2,3,17) + lambda(2,1,11) + lambda(3,3,11) <= 1
 sol(63): lambda(2,3,15) + lambda(3,3,11) + lambda(3,1,12) <= 1
 sol(64): lambda(2,3,17) + lambda(3,3,11) + lambda(3,1,12) <= 1
 sol(65): lambda(2,3,17) + lambda(3,3,16) + lambda(1,1,13) <= 1
 sol(66): lambda(1,1,15) + lambda(3,3,11) + lambda(2,3,12) <= 1
 sol(67): lambda(1,1,15) + lambda(3,2,13) + lambda(2,3,15) <= 1
 sol(68): lambda(3,1,17) + lambda(3,3,11) + lambda(2,3,12) <= 1
 sol(69): lambda(2,3,14) + lambda(3,3,11) + lambda(1,1,13) <= 1
 sol(70): lambda(2,3,15) + lambda(3,3,11) + lambda(1,1,13) <= 1
 sol(71): lambda(2,3,17) + lambda(1,1,15) + lambda(3,2,13) <= 1
 sol(72): lambda(2,3,14) + lambda(3,2,13) + lambda(1,1,13) <= 1
 sol(73): lambda(2,3,17) + lambda(1,1,15) + lambda(2,1,11) + lambda(3,3,11)
   <= 1
 sol(74): lambda(2,3,17) + lambda(3,3,11) + lambda(1,1,13) <= 1
 sol(75): lambda(2,3,17) + lambda(3,3,16) + lambda(1,1,15) + lambda(3,2,13)
   + lambda(2,1,11) <= 1
 sol(76): lambda(2,3,17) + lambda(2,1,13) + lambda(3,3,11) <= 1
 sol(77): lambda(3,1,17) + lambda(2,3,14) + lambda(3,3,11) <= 1
 sol(78): lambda(1,1,15) + lambda(2,3,15) + lambda(3,3,11) <= 1
 sol(79): lambda(2,3,17) + lambda(3,1,17) + lambda(3,3,11) <= 1
 sol(80): lambda(2,3,17) + lambda(3,1,10) + lambda(1,1,13) <= 1
 sol(81): lambda(2,3,17) + lambda(3,3,11) <= 1
 sol(82): lambda(1,1,15) + lambda(2,3,14) + lambda(3,3,11) <= 1
 sol(83): lambda(1,1,15) + lambda(2,3,14) + lambda(3,1,12) <= 1
 sol(84): lambda(1,1,15) + lambda(2,3,15) + lambda(3,1,12) <= 1
 sol(85): lambda(2,3,17) + lambda(2,1,13) <= 1
 sol(86): lambda(2,3,14) + lambda(1,1,13) <= 1
 sol(87): lambda(2,3,17) + lambda(3,1,17) + lambda(3,2,13) <= 1
 sol(88): lambda(3,2,13) + lambda(2,3,15) + lambda(1,1,13) <= 1
 sol(89): lambda(3,1,17) + lambda(3,2,13) + lambda(2,3,12) <= 1
 sol(90): lambda(2,3,17) + lambda(3,2,13) + lambda(1,1,13) <= 1
 sol(91): lambda(2,3,15) + lambda(3,1,12) <= 1
 sol(92): lambda(3,1,17) + lambda(2,3,15) + lambda(3,3,11) <= 1
 sol(93): lambda(3,3,16) + lambda(3,2,13) + lambda(2,3,12) + lambda(1,1,13)
   <= 1
 sol(94): lambda(3,3,16) + lambda(2,3,12) + lambda(1,1,13) <= 1
 sol(95): lambda(3,3,16) + lambda(1,1,15) + lambda(3,2,13) + lambda(2,3,12)
   <= 1
 sol(96): lambda(3,3,16) + lambda(1,1,15) + lambda(2,3,12) <= 1
 sol(97): lambda(2,3,17) + lambda(3,2,13) + lambda(2,1,11) <= 1
 sol(98): lambda(3,1,17) + lambda(3,3,16) + lambda(3,2,13) + lambda(2,3,12)
   <= 1
 sol(99): lambda(2,3,17) + lambda(1,1,15) + lambda(3,1,12) <= 1
 sol(100): lambda(3,1,17) + lambda(3,3,16) + lambda(2,3,12) <= 1
 sol(101): lambda(2,3,17) + lambda(3,1,10) <= 1
 sol(102): lambda(2,3,17) + lambda(3,1,17) + lambda(3,2,13)
   + lambda(2,1,13) <= 1
Bounds
 lambda(2,3,17) <= 1
 lambda(3,1,17) <= 1
 lambda(3,3,16) <= 1
 lambda(1,1,15) <= 1
 lambda(2,3,14) <= 1
 lambda(3,2,13) <= 1
 lambda(2,3,15) <= 1
 lambda(2,1,13) <= 1
 lambda(2,1,11) <= 1
 lambda(3,3,11) <= 1
 lambda(3,1,12) <= 1
 lambda(3,1,10) <= 1
 lambda(2,3,12) <= 1
 lambda(1,1,13) <= 1
End
