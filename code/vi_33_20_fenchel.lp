\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.1028225806451613 lambda(3,2,0) - 0.1431451612903226 lambda(3,2,3)
   - 0.0806451612903226 lambda(2,2,6) - 0.1431451612903226 lambda(1,3,6)
   - 0.0806451612903226 lambda(3,1,9) - 0.1431451612903226 lambda(3,1,7)
   - 0.4798387096774193 lambda(1,3,3) - 0.4173387096774194 lambda(2,3,9)
   - 0.4173387096774194 lambda(2,2,0) - 0.1028225806451613 lambda(3,1,4)
   - 0.4798387096774193 lambda(1,2,0) - 0.2056451612903226 lambda(2,2,8)
   - 0.0806451612903226 lambda(3,2,5) - 0.1028225806451613 lambda(3,3,3)
   - 0.3366935483870968 lambda(1,3,8) - 0.4173387096774194 lambda(2,1,5)
   - 0.1431451612903226 lambda(3,3,6) - 0.0806451612903226 lambda(3,3,8)
   - 0.439516129032258 lambda(1,1,8) - 0.0403225806451613 lambda(1,1,9)
   - 0.0403225806451613 lambda(1,2,6) - 0.0403225806451613 lambda(1,3,9)
   - 0.6733870967741935 lambda(3,2,9) - 0.2963709677419354 lambda(2,1,9)
   - 0.2963709677419354 lambda(2,2,3) - 0.1431451612903226 lambda(1,2,3)
   - 0.3366935483870968 lambda(1,2,5)
Subject To
 sol(0): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,6)
   + lambda(3,1,9) <= 1
 sol(1): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,6)
   + lambda(3,1,7) <= 1
 sol(2): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(2,3,9) <= 1
 sol(3): lambda(3,1,7) + lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0)
   + lambda(3,1,4) <= 1
 sol(4): lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,9) + lambda(1,2,0)
   + lambda(2,2,8) <= 1
 sol(5): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(1,3,8) <= 1
 sol(6): lambda(3,2,0) + lambda(3,1,9) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(7): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,9) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) <= 1
 sol(8): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,8) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(9): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(2,2,8)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(10): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(11): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(12): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(13): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(2,3,9)
   + lambda(3,3,6) + lambda(1,2,6) <= 1
 sol(14): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(15): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,6) + lambda(2,1,5)
   + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(16): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(17): lambda(3,2,0) + lambda(1,3,3) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(18): lambda(3,2,0) + lambda(1,3,3) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(19): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,2,5)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(20): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,7)
   + lambda(2,2,8) <= 1
 sol(21): lambda(3,2,0) + lambda(2,2,6) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(22): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(3,1,4) + lambda(1,3,8) <= 1
 sol(23): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(24): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(2,2,8) + lambda(1,3,9) <= 1
 sol(25): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(26): lambda(3,2,0) + lambda(2,2,6) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(27): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9)
   + lambda(1,3,3) + lambda(1,3,8) <= 1
 sol(28): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,3)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(29): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(30): lambda(3,2,0) + lambda(3,1,7) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(3,1,4) + lambda(2,2,3) <= 1
 sol(31): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,9) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(1,2,3) <= 1
 sol(32): lambda(2,2,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(33): lambda(2,2,6) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(34): lambda(2,2,6) + lambda(3,1,7) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) <= 1
 sol(35): lambda(3,1,7) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(36): lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(37): lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(38): lambda(1,3,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,3) <= 1
 sol(39): lambda(3,1,9) + lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(40): lambda(2,2,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,8) <= 1
 sol(41): lambda(2,2,6) + lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,9) <= 1
 sol(42): lambda(2,2,6) + lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,8) <= 1
 sol(43): lambda(1,3,6) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(44): lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,9) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(3,3,3) <= 1
 sol(45): lambda(3,1,7) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,2,6) <= 1
 sol(46): lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(47): lambda(1,3,3) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(48): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(1,3,8) <= 1
 sol(49): lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(50): lambda(1,3,6) + lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) <= 1
 sol(51): lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(52): lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(53): lambda(2,2,6) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(54): lambda(2,2,6) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(55): lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(1,3,8) <= 1
 sol(56): lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(57): lambda(1,3,6) + lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,3) <= 1
 sol(58): lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(59): lambda(3,2,3) + lambda(1,3,3) + lambda(2,3,9) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(60): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(61): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(1,2,0) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(62): lambda(2,3,9) + lambda(1,2,0) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(63): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9)
   + lambda(1,3,3) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(64): lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,7) + lambda(1,2,0)
   + lambda(2,2,8) <= 1
 sol(65): lambda(2,2,6) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(66): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(67): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7)
   + lambda(3,3,6) + lambda(1,3,9) <= 1
 sol(68): lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(69): lambda(3,2,3) + lambda(1,3,3) + lambda(2,3,9) + lambda(1,2,0)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(70): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(71): lambda(2,2,6) + lambda(3,1,7) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(72): lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(73): lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(74): lambda(2,2,6) + lambda(3,1,7) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(75): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(2,2,8)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(76): lambda(3,1,7) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(77): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(78): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(79): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(80): lambda(3,2,3) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(81): lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,9) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(82): lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(83): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(84): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(85): lambda(3,2,0) + lambda(3,1,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(86): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9)
   + lambda(1,2,5) <= 1
 sol(87): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(88): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(89): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(90): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(91): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(2,3,9)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) <= 1
 sol(92): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(93): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(94): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(95): lambda(2,2,6) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(96): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(1,3,8) <= 1
 sol(97): lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(98): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(99): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(100): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(101): lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(2,1,9) <= 1
 sol(102): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,9)
   + lambda(2,2,8) <= 1
 sol(103): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(104): lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(1,3,8) <= 1
 sol(105): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(106): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(3,1,4) + lambda(1,3,9) <= 1
 sol(107): lambda(3,2,0) + lambda(1,3,3) + lambda(3,2,5) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(108): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(109): lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(1,3,9) <= 1
 sol(110): lambda(3,2,0) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(111): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(112): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(113): lambda(3,1,7) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(114): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(115): lambda(3,2,0) + lambda(1,3,3) + lambda(3,1,4) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(116): lambda(3,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(117): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(118): lambda(3,2,0) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(119): lambda(2,2,6) + lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(120): lambda(3,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(121): lambda(3,2,0) + lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,9)
   + lambda(3,1,4) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(122): lambda(3,1,9) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) <= 1
 sol(123): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(124): lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(125): lambda(3,2,3) + lambda(1,3,3) + lambda(1,2,0) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(126): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(127): lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(128): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(129): lambda(3,1,9) + lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(130): lambda(2,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(131): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(132): lambda(1,3,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) <= 1
 sol(133): lambda(1,3,6) + lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(134): lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(135): lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(136): lambda(3,1,7) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(137): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(138): lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(139): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(140): lambda(3,2,3) + lambda(3,1,9) + lambda(1,2,0) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(141): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(3,1,4) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(142): lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(143): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(144): lambda(2,2,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(145): lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(146): lambda(3,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(147): lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(148): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(149): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(150): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(151): lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9)
   + lambda(1,2,5) <= 1
 sol(152): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(153): lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(154): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(155): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,3) <= 1
 sol(156): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(157): lambda(3,2,0) + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(158): lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,7)
   + lambda(1,2,0) <= 1
 sol(159): lambda(1,3,6) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(160): lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(161): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(3,1,4) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(162): lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(163): lambda(2,2,6) + lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(164): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(165): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(166): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(2,1,9) <= 1
 sol(167): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(168): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(169): lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,5) <= 1
 sol(170): lambda(2,2,6) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(171): lambda(3,2,0) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(172): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,9)
   + lambda(2,1,5) <= 1
 sol(173): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6)
   + lambda(1,2,3) <= 1
 sol(174): lambda(2,2,6) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(175): lambda(3,1,9) + lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(176): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(177): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(178): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,9) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,3) <= 1
 sol(179): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(180): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(181): lambda(3,1,9) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,6) <= 1
 sol(182): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(183): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,3) <= 1
 sol(184): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(185): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(186): lambda(2,2,6) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(187): lambda(2,2,6) + lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(188): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(3,2,9) <= 1
 sol(189): lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9)
   + lambda(1,2,3) <= 1
 sol(190): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(191): lambda(1,3,3) + lambda(2,2,0) + lambda(3,2,5) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(192): lambda(2,3,9) + lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(193): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(194): lambda(1,3,3) + lambda(2,2,0) + lambda(3,2,5) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(195): lambda(3,1,7) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,9) <= 1
 sol(196): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(197): lambda(3,1,7) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,6) <= 1
 sol(198): lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(199): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(200): lambda(3,2,0) + lambda(2,3,9) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(201): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(202): lambda(3,2,3) + lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,6) + lambda(1,3,9) <= 1
 sol(203): lambda(3,2,0) + lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,7)
   + lambda(3,1,4) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(204): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(205): lambda(3,1,9) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(206): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(207): lambda(1,3,6) + lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(1,2,3) <= 1
 sol(208): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(3,2,9) <= 1
 sol(209): lambda(3,1,9) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(210): lambda(2,2,6) + lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(211): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(212): lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(213): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(214): lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(215): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(216): lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9) <= 1
 sol(217): lambda(3,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(218): lambda(1,3,3) + lambda(1,2,0) + lambda(3,2,5) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(219): lambda(3,1,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(220): lambda(3,2,0) + lambda(2,2,6) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(221): lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(3,2,9) <= 1
 sol(222): lambda(3,1,9) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(223): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(224): lambda(3,2,0) + lambda(2,2,8) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(225): lambda(1,3,3) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(226): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(227): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7)
   + lambda(1,3,3) + lambda(1,3,9) <= 1
 sol(228): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(229): lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(230): lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(231): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(232): lambda(3,1,9) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(233): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(234): lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(235): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(1,2,3) <= 1
 sol(236): lambda(3,1,7) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(237): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(238): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(239): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(240): lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(241): lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(242): lambda(3,2,3) + lambda(1,2,0) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(2,1,9) <= 1
 sol(243): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(244): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(245): lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(246): lambda(3,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(3,3,6)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(247): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(248): lambda(1,3,6) + lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(249): lambda(2,3,9) + lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(250): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,2,0)
   + lambda(3,3,6) + lambda(1,3,9) <= 1
 sol(251): lambda(2,2,6) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(252): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(253): lambda(1,3,6) + lambda(3,1,9) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) <= 1
 sol(254): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(255): lambda(2,2,6) + lambda(3,1,7) + lambda(2,2,0) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(256): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,3) <= 1
 sol(257): lambda(1,3,6) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(258): lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(259): lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(260): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(261): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(262): lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(263): lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(264): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(265): lambda(3,2,0) + lambda(3,2,3) + lambda(2,3,9) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(266): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(267): lambda(3,2,3) + lambda(2,3,9) + lambda(1,2,0) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(3,2,9) <= 1
 sol(268): lambda(3,2,0) + lambda(3,1,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9)
   + lambda(1,2,3) <= 1
 sol(269): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(270): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,5) <= 1
 sol(271): lambda(2,2,6) + lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(272): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(2,1,9) <= 1
 sol(273): lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(274): lambda(3,1,7) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(275): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(276): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) <= 1
 sol(277): lambda(1,3,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) <= 1
 sol(278): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(279): lambda(3,2,3) + lambda(1,3,3) + lambda(1,2,0) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(280): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(281): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(282): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(283): lambda(3,2,0) + lambda(3,2,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9)
   + lambda(2,1,9) <= 1
 sol(284): lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(285): lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(286): lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(287): lambda(2,2,6) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(288): lambda(3,2,3) + lambda(3,1,9) + lambda(1,2,0) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) <= 1
 sol(289): lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(290): lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(291): lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(292): lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(1,3,9) <= 1
 sol(293): lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(294): lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(295): lambda(2,2,6) + lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(296): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(297): lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(298): lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,7) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(3,3,3) <= 1
 sol(299): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(300): lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(301): lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(302): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(303): lambda(3,2,0) + lambda(3,1,7) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(304): lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(305): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7)
   + lambda(1,3,3) + lambda(1,3,8) <= 1
 sol(306): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(307): lambda(3,2,0) + lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(308): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(309): lambda(3,2,3) + lambda(3,1,9) + lambda(2,3,9) + lambda(1,2,0)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) <= 1
 sol(310): lambda(2,2,6) + lambda(2,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(311): lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(1,2,0) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(312): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,2,5) <= 1
 sol(313): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(1,3,8) + lambda(2,1,5) <= 1
 sol(314): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(315): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(316): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(317): lambda(3,2,0) + lambda(1,3,3) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(318): lambda(3,2,0) + lambda(2,3,9) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(319): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(320): lambda(3,1,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(321): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,5) <= 1
 sol(322): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(323): lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(324): lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9) + lambda(1,2,0)
   + lambda(3,2,5) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(325): lambda(3,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(326): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(327): lambda(1,3,3) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(328): lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(329): lambda(1,3,6) + lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(330): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(331): lambda(1,3,6) + lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(332): lambda(1,3,6) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(333): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(334): lambda(3,1,7) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(335): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(336): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(337): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(338): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(339): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(3,1,4) + lambda(1,3,8) <= 1
 sol(340): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(341): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(342): lambda(3,2,3) + lambda(2,2,6) + lambda(1,2,0) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(343): lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(344): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(1,2,0) + lambda(1,3,8) <= 1
 sol(345): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(346): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(347): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(348): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,8) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(349): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(350): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(351): lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(352): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(1,2,3) <= 1
 sol(353): lambda(3,1,7) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(354): lambda(1,2,0) + lambda(3,3,3) + lambda(1,3,8) + lambda(3,3,6)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(355): lambda(3,2,0) + lambda(3,3,3) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(356): lambda(3,2,0) + lambda(3,2,3) + lambda(2,3,9) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6)
   + lambda(3,2,9) <= 1
 sol(357): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(358): lambda(2,2,6) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,3) <= 1
 sol(359): lambda(3,1,9) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(360): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,2,0)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(361): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,3)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(362): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(363): lambda(3,2,3) + lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(364): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(365): lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(366): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(367): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(368): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(369): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(370): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(371): lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(1,2,0) <= 1
 sol(372): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(373): lambda(3,2,0) + lambda(1,3,3) + lambda(3,2,5) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(374): lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(375): lambda(3,1,9) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(376): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(377): lambda(1,3,6) + lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(378): lambda(3,2,3) + lambda(1,2,0) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(379): lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(380): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(381): lambda(3,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(382): lambda(2,3,9) + lambda(1,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(383): lambda(3,1,7) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(2,2,3) <= 1
 sol(384): lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(385): lambda(3,2,3) + lambda(2,2,6) + lambda(1,2,0) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(386): lambda(1,2,0) + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(387): lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(388): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(389): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(390): lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(391): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(392): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(393): lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,8)
   + lambda(1,2,6) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(394): lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0)
   + lambda(3,2,5) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(395): lambda(3,2,0) + lambda(3,1,9) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(396): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(397): lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(398): lambda(3,2,0) + lambda(3,1,9) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(399): lambda(2,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(400): lambda(1,3,6) + lambda(3,1,7) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(401): lambda(3,1,9) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,2,5) <= 1
 sol(402): lambda(3,1,7) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(403): lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(404): lambda(3,1,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(405): lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(406): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,7) + lambda(2,2,8)
   + lambda(3,3,6) + lambda(1,3,9) <= 1
 sol(407): lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(408): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(409): lambda(3,1,9) + lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(410): lambda(1,2,0) + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,3,9) + lambda(3,2,9) + lambda(2,1,9)
   + lambda(1,2,3) <= 1
 sol(411): lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(412): lambda(2,3,9) + lambda(2,2,0) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(413): lambda(3,1,9) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(414): lambda(3,1,9) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(415): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,9) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(1,2,3) <= 1
 sol(416): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(417): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(418): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,9) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(419): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(420): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9)
   + lambda(1,2,5) <= 1
 sol(421): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(422): lambda(3,2,0) + lambda(2,3,9) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(423): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(424): lambda(3,2,0) + lambda(1,3,6) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,3) <= 1
 sol(425): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(426): lambda(3,2,3) + lambda(1,3,6) + lambda(1,2,0) + lambda(2,1,5)
   + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(427): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,6) + lambda(1,2,3) <= 1
 sol(428): lambda(3,1,7) + lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(3,3,3) + lambda(1,2,5) <= 1
 sol(429): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(430): lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(1,3,8) <= 1
 sol(431): lambda(3,2,3) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(432): lambda(3,2,0) + lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(433): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(434): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(435): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(436): lambda(3,1,9) + lambda(2,2,0) + lambda(3,3,3) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(437): lambda(3,2,0) + lambda(2,2,6) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,3) <= 1
 sol(438): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(439): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,3) <= 1
 sol(440): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(441): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,6) + lambda(3,2,9)
   + lambda(1,2,3) <= 1
 sol(442): lambda(3,2,0) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(443): lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(444): lambda(3,2,0) + lambda(1,3,3) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(445): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(2,3,9) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(446): lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,8) + lambda(2,2,3) <= 1
 sol(447): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(1,2,0) + lambda(1,3,8) <= 1
 sol(448): lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(449): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3)
   + lambda(2,2,8) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(450): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9)
   + lambda(3,2,5) + lambda(2,1,5) + lambda(3,3,8) <= 1
 sol(451): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,2,6) + lambda(1,3,9) + lambda(1,2,3) <= 1
 sol(452): lambda(3,2,0) + lambda(1,3,3) + lambda(2,3,9) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(453): lambda(3,1,9) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,6) + lambda(3,3,8) + lambda(1,2,6)
   + lambda(1,2,3) <= 1
 sol(454): lambda(2,3,9) + lambda(1,2,0) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(455): lambda(3,2,0) + lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(456): lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,6) + lambda(3,1,9)
   + lambda(1,2,0) <= 1
 sol(457): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,9) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(3,3,3) + lambda(2,2,3) <= 1
 sol(458): lambda(3,2,0) + lambda(1,3,6) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(459): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(460): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(2,2,3) <= 1
 sol(461): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(2,2,3) <= 1
 sol(462): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(2,2,3) <= 1
 sol(463): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(464): lambda(3,2,0) + lambda(2,2,8) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(2,2,3) <= 1
 sol(465): lambda(3,2,0) + lambda(3,1,9) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(2,2,3) <= 1
 sol(466): lambda(3,1,9) + lambda(1,3,3) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(467): lambda(3,2,0) + lambda(3,1,9) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(468): lambda(1,3,3) + lambda(2,3,9) + lambda(1,2,0) + lambda(3,2,5)
   + lambda(2,1,5) + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) <= 1
 sol(469): lambda(3,2,0) + lambda(3,2,3) + lambda(1,3,8) + lambda(2,1,5)
   + lambda(3,3,6) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(470): lambda(3,1,9) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(1,3,8) <= 1
 sol(471): lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) <= 1
 sol(472): lambda(3,2,0) + lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(473): lambda(3,2,0) + lambda(2,3,9) + lambda(3,3,3) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(2,2,3) <= 1
 sol(474): lambda(1,3,3) + lambda(1,2,0) + lambda(3,2,5) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(475): lambda(2,3,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(3,2,9) + lambda(1,2,3) <= 1
 sol(476): lambda(3,2,0) + lambda(3,2,3) + lambda(3,1,9) + lambda(2,2,8)
   + lambda(3,3,6) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(477): lambda(3,2,0) + lambda(3,1,7) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(478): lambda(3,2,0) + lambda(3,2,5) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(479): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,8) + lambda(1,2,5) <= 1
 sol(480): lambda(3,2,3) + lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,8) <= 1
 sol(481): lambda(3,2,0) + lambda(3,1,9) + lambda(3,1,4) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,5) <= 1
 sol(482): lambda(3,2,3) + lambda(1,3,6) + lambda(3,1,9) + lambda(1,2,0)
   + lambda(2,1,5) <= 1
 sol(483): lambda(3,2,0) + lambda(2,2,6) + lambda(3,1,7) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,3,8) + lambda(1,2,3) <= 1
 sol(484): lambda(3,2,3) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,3,6)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(485): lambda(3,2,0) + lambda(2,3,9) + lambda(3,1,4) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(3,2,9) + lambda(1,2,5) <= 1
 sol(486): lambda(3,2,0) + lambda(3,1,4) + lambda(3,3,3) + lambda(1,3,8)
   + lambda(3,2,9) + lambda(2,1,9) + lambda(1,2,5) <= 1
 sol(487): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(488): lambda(3,2,0) + lambda(3,1,4) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,8) + lambda(1,1,9) + lambda(1,2,5) <= 1
 sol(489): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,9) + lambda(1,2,0)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(490): lambda(3,1,7) + lambda(1,3,3) + lambda(3,1,4) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(491): lambda(3,2,0) + lambda(3,1,7) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(492): lambda(3,1,9) + lambda(3,1,4) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) + lambda(2,2,3) <= 1
 sol(493): lambda(3,2,0) + lambda(3,1,9) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(1,3,8) + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(494): lambda(3,2,3) + lambda(3,1,7) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(495): lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,3) + lambda(1,3,8) <= 1
 sol(496): lambda(3,2,0) + lambda(1,3,6) + lambda(3,2,5) + lambda(3,3,3)
   + lambda(2,1,5) + lambda(3,2,9) + lambda(2,1,9) <= 1
 sol(497): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,3,3)
   + lambda(1,2,0) + lambda(1,3,9) <= 1
 sol(498): lambda(3,2,3) + lambda(2,2,6) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(499): lambda(3,2,3) + lambda(3,1,7) + lambda(1,3,3) + lambda(1,2,0)
   + lambda(2,2,8) + lambda(1,3,9) <= 1
 sol(500): lambda(3,1,7) + lambda(2,3,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(1,2,5) <= 1
 sol(501): lambda(2,2,6) + lambda(3,1,9) + lambda(1,3,3) + lambda(3,1,4)
   + lambda(1,2,0) + lambda(1,3,8) <= 1
 sol(502): lambda(2,2,6) + lambda(3,1,9) + lambda(2,2,0) + lambda(3,1,4)
   + lambda(3,3,3) + lambda(3,3,8) + lambda(1,3,9) <= 1
 sol(503): lambda(3,2,3) + lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,9) <= 1
 sol(504): lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8) + lambda(3,2,5)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(505): lambda(3,1,9) + lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,2,5) + lambda(1,3,8) <= 1
 sol(506): lambda(3,2,3) + lambda(3,1,7) + lambda(1,2,0) + lambda(3,3,6)
   + lambda(1,2,6) + lambda(1,3,9) <= 1
 sol(507): lambda(3,2,0) + lambda(3,1,7) + lambda(2,3,9) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,2,5) <= 1
 sol(508): lambda(3,1,9) + lambda(1,3,3) + lambda(2,3,9) + lambda(2,2,0)
   + lambda(3,1,4) + lambda(3,3,8) <= 1
 sol(509): lambda(3,2,3) + lambda(3,1,7) + lambda(2,3,9) + lambda(1,2,0)
   + lambda(3,3,6) + lambda(1,2,6) <= 1
 sol(510): lambda(3,1,7) + lambda(2,2,0) + lambda(2,2,8) + lambda(3,3,3)
   + lambda(3,3,6) + lambda(1,3,9) + lambda(1,2,5) <= 1
 sol(511): lambda(3,2,3) + lambda(2,2,6) + lambda(3,1,7) + lambda(1,2,0)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
 sol(512): lambda(3,2,3) + lambda(1,3,3) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(3,3,8) + lambda(1,1,8) <= 1
 sol(513): lambda(3,2,3) + lambda(3,1,9) + lambda(1,2,0) + lambda(2,2,8)
   + lambda(1,3,8) + lambda(3,3,6) <= 1
Bounds
 lambda(3,2,0) <= 1
 lambda(3,2,3) <= 1
 lambda(2,2,6) <= 1
 lambda(1,3,6) <= 1
 lambda(3,1,9) <= 1
 lambda(3,1,7) <= 1
 lambda(1,3,3) <= 1
 lambda(2,3,9) <= 1
 lambda(2,2,0) <= 1
 lambda(3,1,4) <= 1
 lambda(1,2,0) <= 1
 lambda(2,2,8) <= 1
 lambda(3,2,5) <= 1
 lambda(3,3,3) <= 1
 lambda(1,3,8) <= 1
 lambda(2,1,5) <= 1
 lambda(3,3,6) <= 1
 lambda(3,3,8) <= 1
 lambda(1,1,8) <= 1
 lambda(1,1,9) <= 1
 lambda(1,2,6) <= 1
 lambda(1,3,9) <= 1
 lambda(3,2,9) <= 1
 lambda(2,1,9) <= 1
 lambda(2,2,3) <= 1
 lambda(1,2,3) <= 1
 lambda(1,2,5) <= 1
End
