\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.75 lambda(3,1,35) - 0.3333333333333323 lambda(1,2,32)
   - 0.5 lambda(2,3,30) - 0.1666666666666688 lambda(1,2,29)
   - 0.5 lambda(3,2,28) - 0.5 lambda(3,3,25) - 0.5 lambda(4,3,28)
   - 0.5 lambda(2,2,24) - 0.5 lambda(4,3,20) - 0.5 lambda(4,1,18)
   - 0.25 lambda(3,1,20) - 0.5 lambda(1,2,20)
Subject To
 sol(0): lambda(3,1,35) <= 1
 sol(1): lambda(3,1,35) + lambda(1,2,32) <= 1
 sol(2): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) <= 1
 sol(3): lambda(3,1,35) + lambda(1,2,29) <= 1
 sol(4): lambda(3,1,35) + lambda(3,2,28) <= 1
 sol(5): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) <= 1
 sol(6): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) <= 1
 sol(7): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28) <= 1
 sol(8): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(4,3,28)
   <= 1
 sol(9): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(4,3,28)
   <= 1
 sol(10): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) + lambda(3,3,25)
   <= 1
 sol(11): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) + lambda(4,3,28)
   <= 1
 sol(12): lambda(3,1,35) + lambda(2,2,24) <= 1
 sol(13): lambda(3,1,35) + lambda(2,3,30) + lambda(2,2,24) <= 1
 sol(14): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(2,2,24)
   <= 1
 sol(15): lambda(3,1,35) + lambda(4,3,28) + lambda(2,2,24) <= 1
 sol(16): lambda(3,1,35) + lambda(1,2,29) + lambda(4,3,28) + lambda(2,2,24)
   <= 1
 sol(17): lambda(3,1,35) + lambda(3,3,25) + lambda(2,2,24) <= 1
 sol(18): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(2,2,24)
   <= 1
 sol(19): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25) + lambda(2,2,24)
   <= 1
 sol(20): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(3,3,25)
   + lambda(2,2,24) <= 1
 sol(21): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(3,3,25)
   <= 1
 sol(22): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(3,3,25)
   + lambda(2,2,24) <= 1
 sol(23): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(2,2,24) <= 1
 sol(24): lambda(3,1,35) + lambda(4,3,20) <= 1
 sol(25): lambda(3,1,35) + lambda(1,2,32) + lambda(4,3,20) <= 1
 sol(26): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) <= 1
 sol(27): lambda(3,1,35) + lambda(2,3,30) + lambda(4,3,20) <= 1
 sol(28): lambda(3,1,35) + lambda(1,2,29) + lambda(4,3,20) <= 1
 sol(29): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(4,3,20)
   <= 1
 sol(30): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(2,2,24)
   <= 1
 sol(31): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(4,1,18)
   <= 1
 sol(32): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(3,3,25)
   + lambda(4,1,18) <= 1
 sol(33): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(3,3,25)
   + lambda(2,2,24) + lambda(4,1,18) <= 1
 sol(34): lambda(2,3,30) + lambda(4,3,20) + lambda(3,1,20) + lambda(1,2,20)
   <= 1
 sol(35): lambda(3,1,35) + lambda(1,2,29) + lambda(2,2,24) <= 1
 sol(36): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(3,3,25)
   + lambda(4,1,18) <= 1
 sol(37): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(3,3,25)
   <= 1
 sol(38): lambda(1,2,32) + lambda(2,3,30) + lambda(2,2,24) + lambda(3,1,20)
   <= 1
 sol(39): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) <= 1
 sol(40): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) + lambda(1,2,20)
   <= 1
 sol(41): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28) + lambda(4,1,18)
   <= 1
 sol(42): lambda(3,1,35) + lambda(1,2,32) + lambda(4,3,28) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(43): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,20) <= 1
 sol(44): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(4,1,18) <= 1
 sol(45): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) <= 1
 sol(46): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25) + lambda(4,1,18)
   <= 1
 sol(47): lambda(3,1,35) + lambda(1,2,29) + lambda(4,3,28) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(48): lambda(1,2,32) + lambda(4,3,28) + lambda(3,1,20) <= 1
 sol(49): lambda(3,1,35) + lambda(3,2,28) + lambda(4,1,18) <= 1
 sol(50): lambda(2,3,30) + lambda(1,2,29) + lambda(2,2,24) + lambda(3,1,20)
   <= 1
 sol(51): lambda(3,1,35) + lambda(2,3,30) + lambda(2,2,24) + lambda(4,1,18)
   <= 1
 sol(52): lambda(3,1,35) + lambda(4,1,18) <= 1
 sol(53): lambda(3,1,35) + lambda(1,2,29) + lambda(4,3,28) <= 1
 sol(54): lambda(2,3,30) + lambda(1,2,29) + lambda(4,3,20) + lambda(3,1,20)
   <= 1
 sol(55): lambda(3,1,35) + lambda(2,3,30) <= 1
 sol(56): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(2,2,24) + lambda(4,1,18) <= 1
 sol(57): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(58): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,28) + lambda(4,1,18)
   <= 1
 sol(59): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28) + lambda(2,2,24)
   <= 1
 sol(60): lambda(3,1,35) + lambda(1,2,32) + lambda(4,3,28) + lambda(4,1,18)
   <= 1
 sol(61): lambda(3,1,35) + lambda(3,3,25) <= 1
 sol(62): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(2,2,24)
   <= 1
 sol(63): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25) + lambda(4,1,18)
   + lambda(1,2,20) <= 1
 sol(64): lambda(3,1,35) + lambda(1,2,32) + lambda(2,2,24) + lambda(4,3,20)
   <= 1
 sol(65): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,28) + lambda(4,1,18)
   + lambda(1,2,20) <= 1
 sol(66): lambda(3,1,35) + lambda(4,3,28) + lambda(2,2,24) + lambda(4,1,18)
   <= 1
 sol(67): lambda(1,2,32) + lambda(3,1,20) <= 1
 sol(68): lambda(3,1,35) + lambda(3,3,25) + lambda(4,1,18) + lambda(1,2,20)
   <= 1
 sol(69): lambda(3,1,35) + lambda(2,3,30) + lambda(4,1,18) + lambda(1,2,20)
   <= 1
 sol(70): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28) + lambda(4,1,18)
   + lambda(1,2,20) <= 1
 sol(71): lambda(3,1,35) + lambda(1,2,20) <= 1
 sol(72): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25) <= 1
 sol(73): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) + lambda(4,1,18)
   + lambda(1,2,20) <= 1
 sol(74): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(4,1,18)
   <= 1
 sol(75): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) + lambda(3,3,25)
   + lambda(4,1,18) + lambda(1,2,20) <= 1
 sol(76): lambda(3,1,35) + lambda(4,3,28) <= 1
 sol(77): lambda(4,3,28) + lambda(3,1,20) + lambda(1,2,20) <= 1
 sol(78): lambda(3,1,35) + lambda(2,2,24) + lambda(4,1,18) <= 1
 sol(79): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(4,1,18) <= 1
 sol(80): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25) + lambda(4,1,18)
   <= 1
 sol(81): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(82): lambda(3,1,35) + lambda(3,2,28) + lambda(4,1,18) + lambda(1,2,20)
   <= 1
 sol(83): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(84): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) + lambda(3,3,25)
   + lambda(1,2,20) <= 1
 sol(85): lambda(1,2,32) + lambda(2,3,30) + lambda(2,2,24) + lambda(4,3,20)
   + lambda(3,1,20) <= 1
 sol(86): lambda(3,1,35) + lambda(1,2,32) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(2,2,24) + lambda(4,1,18) <= 1
 sol(87): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(88): lambda(3,1,35) + lambda(1,2,29) + lambda(2,2,24) + lambda(4,1,18)
   <= 1
 sol(89): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29) + lambda(3,3,25)
   + lambda(2,2,24) + lambda(4,1,18) <= 1
 sol(90): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) + lambda(4,1,18)
   <= 1
 sol(91): lambda(3,1,35) + lambda(4,3,28) + lambda(4,1,18) <= 1
 sol(92): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(4,1,18) <= 1
 sol(93): lambda(1,2,32) + lambda(2,3,30) + lambda(4,3,20) + lambda(3,1,20)
   <= 1
 sol(94): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(95): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28) + lambda(3,3,25)
   + lambda(4,1,18) <= 1
 sol(96): lambda(3,1,35) + lambda(4,1,18) + lambda(1,2,20) <= 1
 sol(97): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(98): lambda(3,1,35) + lambda(3,3,25) + lambda(2,2,24) + lambda(4,1,18)
   <= 1
 sol(99): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(2,2,24) <= 1
 sol(100): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25)
   + lambda(4,3,28) + lambda(1,2,20) <= 1
 sol(101): lambda(3,1,35) + lambda(2,3,30) + lambda(2,2,24)
   + lambda(4,3,20) <= 1
 sol(102): lambda(3,1,35) + lambda(2,2,24) + lambda(4,3,20) <= 1
 sol(103): lambda(3,1,35) + lambda(1,2,29) + lambda(4,1,18) <= 1
 sol(104): lambda(3,1,35) + lambda(1,2,29) + lambda(2,2,24)
   + lambda(4,3,20) <= 1
 sol(105): lambda(3,1,35) + lambda(2,3,30) + lambda(4,1,18) <= 1
 sol(106): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30)
   + lambda(4,1,18) <= 1
 sol(107): lambda(3,1,35) + lambda(1,2,29) + lambda(4,3,28)
   + lambda(4,1,18) <= 1
 sol(108): lambda(3,1,35) + lambda(3,2,28) + lambda(3,3,25)
   + lambda(4,3,28) + lambda(4,1,18) + lambda(1,2,20) <= 1
 sol(109): lambda(3,1,35) + lambda(1,2,32) + lambda(4,1,18) <= 1
 sol(110): lambda(2,3,30) + lambda(3,1,20) + lambda(1,2,20) <= 1
 sol(111): lambda(3,1,35) + lambda(1,2,32) + lambda(4,3,28) <= 1
 sol(112): lambda(1,2,29) + lambda(4,3,28) + lambda(2,2,24)
   + lambda(3,1,20) <= 1
 sol(113): lambda(3,1,35) + lambda(3,2,28) + lambda(1,2,20) <= 1
 sol(114): lambda(2,3,30) + lambda(1,2,29) + lambda(2,2,24)
   + lambda(4,3,20) + lambda(3,1,20) <= 1
 sol(115): lambda(1,2,32) + lambda(2,2,24) + lambda(4,3,20)
   + lambda(3,1,20) <= 1
 sol(116): lambda(3,1,35) + lambda(3,3,25) + lambda(4,1,18) <= 1
 sol(117): lambda(1,2,32) + lambda(4,3,20) + lambda(3,1,20) <= 1
 sol(118): lambda(1,2,32) + lambda(4,3,28) + lambda(2,2,24)
   + lambda(3,1,20) <= 1
 sol(119): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,28) <= 1
 sol(120): lambda(1,2,29) + lambda(4,3,28) + lambda(3,1,20) <= 1
 sol(121): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28)
   + lambda(4,1,18) + lambda(1,2,20) <= 1
 sol(122): lambda(2,3,30) + lambda(1,2,29) + lambda(3,1,20) <= 1
 sol(123): lambda(1,2,32) + lambda(2,3,30) + lambda(3,1,20) <= 1
 sol(124): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28)
   + lambda(4,3,20) + lambda(1,2,20) <= 1
 sol(125): lambda(1,2,32) + lambda(2,2,24) + lambda(3,1,20) <= 1
 sol(126): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,20)
   + lambda(1,2,20) <= 1
 sol(127): lambda(3,1,35) + lambda(4,3,28) + lambda(4,1,18)
   + lambda(1,2,20) <= 1
 sol(128): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28)
   + lambda(1,2,20) <= 1
 sol(129): lambda(3,1,35) + lambda(2,3,30) + lambda(4,3,20)
   + lambda(1,2,20) <= 1
 sol(130): lambda(3,1,35) + lambda(4,3,20) + lambda(1,2,20) <= 1
 sol(131): lambda(3,1,35) + lambda(1,2,29) + lambda(3,3,25)
   + lambda(4,1,18) <= 1
 sol(132): lambda(3,1,35) + lambda(3,2,28) + lambda(4,3,28)
   + lambda(1,2,20) <= 1
 sol(133): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30)
   + lambda(4,3,20) <= 1
 sol(134): lambda(3,1,35) + lambda(2,3,30) + lambda(3,2,28)
   + lambda(4,3,20) <= 1
 sol(135): lambda(3,1,35) + lambda(3,3,25) + lambda(4,3,28)
   + lambda(1,2,20) <= 1
 sol(136): lambda(3,1,35) + lambda(3,3,25) + lambda(1,2,20) <= 1
 sol(137): lambda(3,1,35) + lambda(1,2,32) + lambda(2,2,24) <= 1
 sol(138): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,29)
   + lambda(2,2,24) + lambda(4,3,20) <= 1
 sol(139): lambda(3,1,35) + lambda(2,3,30) + lambda(3,3,25)
   + lambda(1,2,20) <= 1
 sol(140): lambda(3,1,35) + lambda(4,3,28) + lambda(1,2,20) <= 1
 sol(141): lambda(3,1,35) + lambda(2,3,30) + lambda(1,2,20) <= 1
 sol(142): lambda(3,1,35) + lambda(1,2,32) + lambda(4,3,28)
   + lambda(2,2,24) <= 1
 sol(143): lambda(3,1,35) + lambda(1,2,32) + lambda(2,2,24)
   + lambda(4,1,18) <= 1
 sol(144): lambda(3,1,35) + lambda(1,2,32) + lambda(2,3,30)
   + lambda(2,2,24) + lambda(4,3,20) <= 1
Bounds
 lambda(3,1,35) <= 1
 lambda(1,2,32) <= 1
 lambda(2,3,30) <= 1
 lambda(1,2,29) <= 1
 lambda(3,2,28) <= 1
 lambda(3,3,25) <= 1
 lambda(4,3,28) <= 1
 lambda(2,2,24) <= 1
 lambda(4,3,20) <= 1
 lambda(4,1,18) <= 1
 lambda(3,1,20) <= 1
 lambda(1,2,20) <= 1
End
