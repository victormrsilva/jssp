\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.0408163265306123 lambda(3,4,9) - 0.2408163265306126 lambda(2,1,12)
   - 0.3591836734693878 lambda(2,3,13) - 0.044897959183674 lambda(2,1,13)
   - 0.3795918367346939 lambda(1,3,12) - 0.1795918367346939 lambda(3,4,8)
   - 0.2 lambda(1,3,11) - 0.3591836734693878 lambda(2,2,8)
   - 0.0612244897959183 lambda(1,3,13) - 0.4 lambda(4,1,12)
   - 0.1795918367346939 lambda(4,1,11) - 0.0612244897959183 lambda(4,1,10)
   - 0.1183673469387752 lambda(1,1,6) - 0.4 lambda(4,4,9)
   - 0.289795918367347 lambda(2,4,13) - 0.1387755102040817 lambda(3,2,13)
   - 0.2 lambda(3,4,13) - 0.1795918367346939 lambda(4,4,8)
   - 0.4408163265306123 lambda(3,2,12) - 0.2408163265306126 lambda(4,4,12)
   - 0.0367346938775508 lambda(4,4,13) - 0.1387755102040817 lambda(3,3,10)
   - 0.4408163265306123 lambda(3,3,9) - 0.0612244897959183 lambda(1,4,12)
   - 0.0612244897959183 lambda(4,4,7) - 0.2 lambda(1,4,10)
   - 0.0122448979591836 lambda(1,4,14) - 0.3795918367346939 lambda(1,4,11)
Subject To
 sol(0): lambda(3,4,9) + lambda(2,1,12) + lambda(2,3,13) <= 1
 sol(1): lambda(3,4,9) + lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   <= 1
 sol(2): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(3,4,8)
   <= 1
 sol(3): lambda(2,1,12) + lambda(2,1,13) + lambda(3,4,8) + lambda(1,3,11)
   <= 1
 sol(4): lambda(3,4,9) + lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13)
   <= 1
 sol(5): lambda(3,4,9) + lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8)
   <= 1
 sol(6): lambda(3,4,9) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(7): lambda(3,4,9) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11)
   <= 1
 sol(8): lambda(3,4,9) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(9): lambda(3,4,9) + lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(10): lambda(3,4,9) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(11): lambda(3,4,8) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12)
   <= 1
 sol(12): lambda(2,1,13) + lambda(1,3,12) + lambda(3,4,8) + lambda(2,2,8)
   <= 1
 sol(13): lambda(2,1,13) + lambda(3,4,8) + lambda(1,3,11) + lambda(2,2,8)
   <= 1
 sol(14): lambda(3,4,8) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(15): lambda(3,4,9) + lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6)
   <= 1
 sol(16): lambda(2,3,13) + lambda(3,4,8) + lambda(2,2,8) + lambda(1,1,6)
   <= 1
 sol(17): lambda(2,1,12) + lambda(1,3,11) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) <= 1
 sol(18): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,9)
   + lambda(3,4,13) <= 1
 sol(19): lambda(1,3,11) + lambda(4,1,12) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) <= 1
 sol(20): lambda(2,1,12) + lambda(1,3,11) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(4,4,8) <= 1
 sol(21): lambda(2,3,13) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(22): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(3,4,13)
   + lambda(4,4,8) <= 1
 sol(23): lambda(1,3,11) + lambda(4,1,11) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) <= 1
 sol(24): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,4,13) <= 1
 sol(25): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) <= 1
 sol(26): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) <= 1
 sol(27): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) <= 1
 sol(28): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,13) <= 1
 sol(29): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,10) <= 1
 sol(30): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(31): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(32): lambda(1,3,12) + lambda(4,1,12) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) <= 1
 sol(33): lambda(2,1,12) + lambda(1,3,12) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) <= 1
 sol(34): lambda(3,4,8) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10)
   <= 1
 sol(35): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) <= 1
 sol(36): lambda(3,4,9) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12)
   <= 1
 sol(37): lambda(2,1,12) + lambda(1,3,12) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) <= 1
 sol(38): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) <= 1
 sol(39): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) <= 1
 sol(40): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,12) <= 1
 sol(41): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) <= 1
 sol(42): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) <= 1
 sol(43): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) <= 1
 sol(44): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(2,4,13) + lambda(3,2,13) <= 1
 sol(45): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) <= 1
 sol(46): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) <= 1
 sol(47): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) <= 1
 sol(48): lambda(1,1,6) + lambda(2,4,13) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(4,4,7) <= 1
 sol(49): lambda(1,3,11) + lambda(4,1,12) + lambda(2,4,13) + lambda(3,2,12)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(50): lambda(1,3,12) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(51): lambda(2,1,12) + lambda(1,3,11) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(52): lambda(1,3,11) + lambda(4,1,11) + lambda(2,4,13) + lambda(3,2,12)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(53): lambda(1,3,11) + lambda(4,1,10) + lambda(2,4,13) + lambda(3,2,12)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(54): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,4,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(55): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(56): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(57): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(58): lambda(2,3,13) + lambda(1,1,6) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(59): lambda(2,3,13) + lambda(1,1,6) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(60): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(3,4,13)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(61): lambda(2,3,13) + lambda(1,1,6) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(62): lambda(1,1,6) + lambda(2,4,13) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(3,3,9) <= 1
 sol(63): lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) <= 1
 sol(64): lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) <= 1
 sol(65): lambda(3,4,9) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10)
   <= 1
 sol(66): lambda(1,1,6) + lambda(4,4,9) + lambda(2,4,13) + lambda(3,2,12)
   + lambda(3,3,9) <= 1
 sol(67): lambda(2,2,8) + lambda(1,1,6) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,9) <= 1
 sol(68): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) <= 1
 sol(69): lambda(2,3,13) + lambda(1,1,6) + lambda(3,2,12) + lambda(4,4,13)
   + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(70): lambda(1,3,12) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(71): lambda(2,2,8) + lambda(1,1,6) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(4,4,7) <= 1
 sol(72): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(4,4,7) <= 1
 sol(73): lambda(3,4,9) + lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(74): lambda(2,2,8) + lambda(1,1,6) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(75): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(76): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,10) <= 1
 sol(77): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(78): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(79): lambda(2,1,12) + lambda(1,3,13) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(80): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(81): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) <= 1
 sol(82): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(83): lambda(1,3,13) + lambda(4,1,12) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(84): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(85): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,12) <= 1
 sol(86): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(87): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(88): lambda(2,1,12) + lambda(2,1,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(89): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(90): lambda(2,1,12) + lambda(1,3,11) + lambda(2,4,13) + lambda(3,2,12)
   + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(91): lambda(2,1,12) + lambda(2,1,13) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(92): lambda(2,1,12) + lambda(1,3,13) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(93): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(94): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(95): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(96): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(97): lambda(2,1,12) + lambda(2,3,13) + lambda(3,4,13) + lambda(4,4,8)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(98): lambda(2,1,12) + lambda(2,1,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(99): lambda(1,3,12) + lambda(4,1,12) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(100): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) <= 1
 sol(101): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(102): lambda(1,3,13) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(103): lambda(2,1,13) + lambda(2,2,8) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(104): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,10)
   + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(105): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(106): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(107): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(108): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(109): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,13) <= 1
 sol(110): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(111): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(112): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(113): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(114): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,12)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(115): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(116): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(117): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(118): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(119): lambda(3,4,9) + lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   <= 1
 sol(120): lambda(2,1,12) + lambda(2,1,13) + lambda(3,4,8) + lambda(1,3,13)
   <= 1
 sol(121): lambda(2,1,12) + lambda(2,3,13) + lambda(3,4,13) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(122): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11)
   <= 1
 sol(123): lambda(3,4,8) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11)
   <= 1
 sol(124): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(125): lambda(2,1,12) + lambda(2,1,13) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(126): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(127): lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(128): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(129): lambda(2,1,12) + lambda(2,1,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(130): lambda(2,3,13) + lambda(4,1,12) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(131): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(132): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(133): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(134): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(135): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(2,4,13) + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12)
   <= 1
 sol(136): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(137): lambda(2,1,12) + lambda(1,3,11) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) <= 1
 sol(138): lambda(2,1,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11)
   <= 1
 sol(139): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(140): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(141): lambda(2,3,13) + lambda(4,1,10) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(142): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(143): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(144): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(145): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(146): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,12)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(147): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(148): lambda(2,1,12) + lambda(1,3,13) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(149): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(150): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(151): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(152): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,4,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(153): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,12)
   + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(154): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(155): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(156): lambda(2,3,13) + lambda(4,1,12) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(157): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(158): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(159): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(160): lambda(1,3,12) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(161): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(162): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(163): lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(164): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(165): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(166): lambda(2,3,13) + lambda(4,1,12) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(167): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(168): lambda(1,3,12) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(169): lambda(3,4,9) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(170): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(171): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(172): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(173): lambda(2,1,12) + lambda(2,3,13) + lambda(3,4,13) + lambda(4,4,8)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(174): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,10) + lambda(1,4,11) <= 1
 sol(175): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(176): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(177): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(178): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(179): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(180): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(181): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(182): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(183): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(184): lambda(2,3,13) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(185): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(186): lambda(2,1,12) + lambda(1,3,12) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(187): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(188): lambda(2,3,13) + lambda(1,1,6) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(189): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(190): lambda(3,4,8) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(191): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,4,13)
   + lambda(4,4,8) <= 1
 sol(192): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(193): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(194): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,10)
   + lambda(1,4,11) <= 1
 sol(195): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(196): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(197): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(198): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(199): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(200): lambda(2,1,12) + lambda(2,3,13) + lambda(3,4,8) <= 1
 sol(201): lambda(3,4,9) + lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(202): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(203): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(204): lambda(3,4,9) + lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8)
   <= 1
 sol(205): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(206): lambda(3,4,9) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(207): lambda(1,3,12) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,11)
   <= 1
 sol(208): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(209): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(210): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(211): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(212): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(213): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(214): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(215): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,10) <= 1
 sol(216): lambda(2,1,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(3,3,10)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(217): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,10) <= 1
 sol(218): lambda(2,1,12) + lambda(2,1,13) + lambda(3,2,13) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(219): lambda(2,1,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(220): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(221): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(222): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(223): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,12)
   + lambda(3,3,10) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(224): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(225): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(226): lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(227): lambda(2,1,12) + lambda(2,3,13) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(228): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(229): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,10) <= 1
 sol(230): lambda(2,1,12) + lambda(1,3,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(231): lambda(2,1,13) + lambda(3,4,8) + lambda(2,2,8) + lambda(1,3,13)
   <= 1
 sol(232): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) <= 1
 sol(233): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(234): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(235): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,9)
   + lambda(1,4,11) <= 1
 sol(236): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(237): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(238): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(239): lambda(2,1,12) + lambda(1,3,12) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(240): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(241): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(242): lambda(1,3,12) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(243): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(244): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,10) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(245): lambda(2,3,13) + lambda(4,1,11) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(246): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(247): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(248): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(249): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(250): lambda(2,1,12) + lambda(2,1,13) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(251): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(252): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(253): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(254): lambda(2,1,12) + lambda(2,1,13) + lambda(3,2,12) + lambda(3,3,9)
   + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,14)
   + lambda(1,4,11) <= 1
 sol(255): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(256): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14)
   <= 1
 sol(257): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(258): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12)
   <= 1
 sol(259): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(260): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(261): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) <= 1
 sol(262): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,10) <= 1
 sol(263): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) <= 1
 sol(264): lambda(2,3,13) + lambda(4,1,11) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(265): lambda(2,1,13) + lambda(2,2,8) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) <= 1
 sol(266): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(267): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(4,4,7) <= 1
 sol(268): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(269): lambda(3,4,8) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(270): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,9)
   + lambda(1,4,11) <= 1
 sol(271): lambda(2,1,12) + lambda(1,3,11) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) <= 1
 sol(272): lambda(2,1,12) + lambda(1,3,13) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(273): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,12) + lambda(3,3,10) + lambda(1,4,11) <= 1
 sol(274): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,8) + lambda(3,2,12)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(275): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,4,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(276): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(277): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(278): lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(279): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(280): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,9)
   + lambda(1,4,11) <= 1
 sol(281): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,10) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(282): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11)
   <= 1
 sol(283): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11)
   <= 1
 sol(284): lambda(2,3,13) + lambda(3,4,8) + lambda(2,2,8) + lambda(4,1,10)
   <= 1
 sol(285): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(286): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,11) <= 1
 sol(287): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(288): lambda(2,3,13) + lambda(1,1,6) + lambda(3,2,12) + lambda(4,4,12)
   + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(289): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(290): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(291): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12)
   <= 1
 sol(292): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7) <= 1
 sol(293): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,2,13)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(294): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,10) + lambda(4,4,7) <= 1
 sol(295): lambda(1,3,12) + lambda(4,1,11) + lambda(2,4,13) + lambda(4,4,8)
   + lambda(3,2,12) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(296): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,2,12) + lambda(4,4,12) + lambda(3,3,9) + lambda(4,4,7)
   + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(297): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(298): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,12) + lambda(3,3,10) <= 1
 sol(299): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,10) <= 1
 sol(300): lambda(1,3,11) + lambda(4,1,12) + lambda(4,4,9) + lambda(2,4,13)
   + lambda(3,2,12) <= 1
 sol(301): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(302): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,13) + lambda(3,3,9) <= 1
 sol(303): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(2,4,13) + lambda(3,2,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(304): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(305): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,13)
   + lambda(3,4,13) + lambda(4,4,8) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(306): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(307): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,10) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(308): lambda(2,1,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,9) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11) <= 1
 sol(309): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(310): lambda(2,1,12) + lambda(2,3,13) + lambda(4,4,9) + lambda(3,2,12)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(311): lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(312): lambda(2,3,13) + lambda(2,2,8) + lambda(1,1,6) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,10) <= 1
 sol(313): lambda(2,1,13) + lambda(1,3,11) + lambda(2,2,8) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(314): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,11) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(315): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(316): lambda(3,4,9) + lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11)
   <= 1
 sol(317): lambda(2,2,8) + lambda(1,1,6) + lambda(2,4,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(3,3,10) <= 1
 sol(318): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,4,13)
   + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(319): lambda(2,3,13) + lambda(2,2,8) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(4,4,12) + lambda(3,3,10) + lambda(1,4,11) <= 1
 sol(320): lambda(3,4,9) + lambda(1,3,11) + lambda(2,2,8) + lambda(4,1,12)
   <= 1
 sol(321): lambda(2,1,12) + lambda(2,1,13) + lambda(3,2,13) + lambda(4,4,8)
   + lambda(3,3,10) + lambda(1,4,12) + lambda(1,4,14) + lambda(1,4,11)
   <= 1
 sol(322): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(2,4,13) + lambda(3,2,13) + lambda(3,3,9) <= 1
 sol(323): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(324): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(325): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(3,2,13)
   + lambda(4,4,8) + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(326): lambda(2,3,13) + lambda(2,2,8) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,9) + lambda(1,4,12) <= 1
 sol(327): lambda(2,1,13) + lambda(2,2,8) + lambda(1,3,13) + lambda(4,4,9)
   + lambda(3,2,13) + lambda(4,4,12) + lambda(3,3,9) <= 1
 sol(328): lambda(2,1,13) + lambda(1,3,12) + lambda(2,2,8) + lambda(3,4,13)
   + lambda(3,3,9) + lambda(4,4,7) + lambda(1,4,10) + lambda(1,4,11) <= 1
 sol(329): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(330): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,11)
   + lambda(3,2,12) + lambda(4,4,13) + lambda(4,4,7) + lambda(1,4,10) <= 1
 sol(331): lambda(1,3,12) + lambda(2,2,8) + lambda(4,1,11) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(332): lambda(2,3,13) + lambda(2,2,8) + lambda(4,4,9) + lambda(3,2,13)
   + lambda(4,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(333): lambda(2,3,13) + lambda(2,2,8) + lambda(4,1,12) + lambda(4,4,9)
   + lambda(3,4,13) + lambda(3,3,10) + lambda(1,4,12) <= 1
 sol(334): lambda(2,1,12) + lambda(2,1,13) + lambda(1,3,12)
   + lambda(3,4,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
 sol(335): lambda(2,2,8) + lambda(1,3,13) + lambda(4,1,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,12)
   + lambda(1,4,11) <= 1
 sol(336): lambda(2,1,12) + lambda(2,3,13) + lambda(3,2,13)
   + lambda(4,4,12) + lambda(3,3,10) + lambda(4,4,7) + lambda(1,4,10)
   + lambda(1,4,11) <= 1
 sol(337): lambda(2,1,12) + lambda(1,3,12) + lambda(2,4,13)
   + lambda(3,2,13) + lambda(4,4,8) + lambda(3,3,9) + lambda(1,4,11) <= 1
Bounds
 lambda(3,4,9) <= 1
 lambda(2,1,12) <= 1
 lambda(2,3,13) <= 1
 lambda(2,1,13) <= 1
 lambda(1,3,12) <= 1
 lambda(3,4,8) <= 1
 lambda(1,3,11) <= 1
 lambda(2,2,8) <= 1
 lambda(1,3,13) <= 1
 lambda(4,1,12) <= 1
 lambda(4,1,11) <= 1
 lambda(4,1,10) <= 1
 lambda(1,1,6) <= 1
 lambda(4,4,9) <= 1
 lambda(2,4,13) <= 1
 lambda(3,2,13) <= 1
 lambda(3,4,13) <= 1
 lambda(4,4,8) <= 1
 lambda(3,2,12) <= 1
 lambda(4,4,12) <= 1
 lambda(4,4,13) <= 1
 lambda(3,3,10) <= 1
 lambda(3,3,9) <= 1
 lambda(1,4,12) <= 1
 lambda(4,4,7) <= 1
 lambda(1,4,10) <= 1
 lambda(1,4,14) <= 1
 lambda(1,4,11) <= 1
End
