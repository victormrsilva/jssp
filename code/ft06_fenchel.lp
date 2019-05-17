\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.2784810126582278 lambda(2,2,8) - 0.0379746835443056 lambda(3,3,11)
   - 0.5443037974683533 lambda(1,1,8) - 0.0126582278481015 lambda(6,4,11)
   - 0.4430379746835444 lambda(6,1,15) - 0.0759493670886049 lambda(3,3,13)
   - 0.4050632911392399 lambda(4,1,13) - 0.075949367088599 lambda(6,6,9)
   - 0.1518987341772155 lambda(1,1,14) - 0.037974683544312 lambda(6,6,8)
   - 0.1518987341772155 lambda(1,3,8) - 0.1518987341772155 lambda(1,1,9)
   - 0.0126582278481015 lambda(6,2,8) - 0.0126582278481018 lambda(6,2,13)
   - 0.4177215189873413 lambda(2,5,13) - 0.2784810126582276 lambda(2,3,11)
   - 0.1518987341772153 lambda(3,4,14) - 0.1139240506329035 lambda(3,6,9)
   - 0.4050632911392399 lambda(4,2,8) - 0.139240506329114 lambda(3,6,11)
   - 0.1772151898734255 lambda(3,6,10) - 0.0253164556962039 lambda(5,3,13)
   - 0.3037974683544311 lambda(3,6,15) - 0.0126582278481015 lambda(6,6,14)
   - 0.2911392405063291 lambda(6,4,14) - 0.0126582278481013 lambda(5,3,15)
   - 0.0506329113924072 lambda(3,4,11) - 0.3037974683544312 lambda(4,3,14)
   - 0.4177215189873413 lambda(2,3,8) - 0.0253164556962033 lambda(2,3,13)
   - 0.2784810126582279 lambda(4,1,8) - 0.2658227848101258 lambda(5,2,13)
   - 0.1518987341772155 lambda(4,2,13) - 0.2911392405063291 lambda(1,2,15)
   - 0.1518987341772155 lambda(1,3,10) - 0.2911392405063291 lambda(6,2,11)
   - 0.1518987341772153 lambda(3,3,9)
Subject To
 sol(0): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(6,1,15) <= 1
 sol(1): lambda(2,2,8) + lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,6,9) <= 1
 sol(2): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13)
   + lambda(6,6,9) <= 1
 sol(3): lambda(2,2,8) + lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,6,8) <= 1
 sol(4): lambda(2,2,8) + lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,6,8) <= 1
 sol(5): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13)
   + lambda(6,6,8) <= 1
 sol(6): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(7): lambda(2,2,8) + lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(8): lambda(2,2,8) + lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(9): lambda(2,2,8) + lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(10): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(11): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) <= 1
 sol(12): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) <= 1
 sol(13): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) <= 1
 sol(14): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) <= 1
 sol(15): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(16): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) <= 1
 sol(17): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) <= 1
 sol(18): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) <= 1
 sol(19): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) <= 1
 sol(20): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) <= 1
 sol(21): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) <= 1
 sol(22): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(5,3,13) + lambda(3,6,15) <= 1
 sol(23): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(6,6,14) <= 1
 sol(24): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(6,4,14) <= 1
 sol(25): lambda(2,2,8) + lambda(1,1,8) + lambda(6,1,15) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(26): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(27): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(28): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14)
   + lambda(6,6,8) <= 1
 sol(29): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(30): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(3,4,11) <= 1
 sol(31): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(1,1,14) + lambda(6,6,14) <= 1
 sol(32): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(6,4,14) <= 1
 sol(33): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,10) + lambda(5,3,15) <= 1
 sol(34): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,10) + lambda(4,3,14) <= 1
 sol(35): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,10) + lambda(5,3,13) <= 1
 sol(36): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(37): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(38): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(39): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(40): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(41): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,9) + lambda(5,3,15) <= 1
 sol(42): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,9) + lambda(5,3,13) <= 1
 sol(43): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(44): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(45): lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(2,3,8) <= 1
 sol(46): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(47): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(48): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(49): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(50): lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) <= 1
 sol(51): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(52): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(53): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(54): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(55): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(56): lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(57): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(58): lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(59): lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(60): lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(61): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(62): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(63): lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(64): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(65): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(66): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(67): lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(68): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(69): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,6,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(70): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,6,10)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(71): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(72): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(73): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,11) + lambda(2,3,13) <= 1
 sol(74): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(75): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,4,14) <= 1
 sol(76): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,11) <= 1
 sol(77): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,10) + lambda(2,3,13) <= 1
 sol(78): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,10) <= 1
 sol(79): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,9) <= 1
 sol(80): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,9) <= 1
 sol(81): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(82): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(83): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(84): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(85): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(86): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(1,3,10) <= 1
 sol(87): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(1,3,10) <= 1
 sol(88): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(89): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(90): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(91): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) <= 1
 sol(92): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) <= 1
 sol(93): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) <= 1
 sol(94): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) <= 1
 sol(95): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(2,3,13) <= 1
 sol(96): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,9) + lambda(4,3,14) <= 1
 sol(97): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) <= 1
 sol(98): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) <= 1
 sol(99): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(2,3,13) <= 1
 sol(100): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(2,3,13) <= 1
 sol(101): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) <= 1
 sol(102): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,4,14) <= 1
 sol(103): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,6,14) <= 1
 sol(104): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(5,3,15) <= 1
 sol(105): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(4,3,14) <= 1
 sol(106): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(107): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(108): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(109): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(110): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,4,14) <= 1
 sol(111): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,6,14) <= 1
 sol(112): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,4,14) <= 1
 sol(113): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(5,3,15) <= 1
 sol(114): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(4,3,14) <= 1
 sol(115): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(5,3,13) <= 1
 sol(116): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(117): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(118): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(119): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(120): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(121): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(122): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(123): lambda(2,2,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(124): lambda(2,2,8) + lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(125): lambda(2,2,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(126): lambda(2,2,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(127): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(128): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(1,2,15) <= 1
 sol(129): lambda(3,3,11) + lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(130): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(1,2,15) <= 1
 sol(131): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(4,2,13) <= 1
 sol(132): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(133): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(134): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(1,2,15) <= 1
 sol(135): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(5,2,13) <= 1
 sol(136): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(137): lambda(2,2,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(138): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(139): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(140): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(141): lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(142): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(6,2,11) <= 1
 sol(143): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(5,2,13) <= 1
 sol(144): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,2,13) <= 1
 sol(145): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(146): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(5,2,13) <= 1
 sol(147): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(148): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(1,2,15) <= 1
 sol(149): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(4,2,13) <= 1
 sol(150): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(1,2,15) <= 1
 sol(151): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,2,13) <= 1
 sol(152): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(153): lambda(1,1,8) + lambda(6,1,15) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(154): lambda(1,1,8) + lambda(1,1,14) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(155): lambda(1,1,8) + lambda(1,1,14) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(156): lambda(1,1,8) + lambda(1,1,14) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(157): lambda(1,1,8) + lambda(4,1,13) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(158): lambda(1,1,8) + lambda(4,1,13) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(159): lambda(1,1,8) + lambda(4,1,13) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(160): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(161): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(162): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(163): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(164): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(165): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(166): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(167): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(168): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(169): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(170): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(171): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(172): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) <= 1
 sol(173): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(174): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(175): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(176): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(177): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) <= 1
 sol(178): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) <= 1
 sol(179): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) <= 1
 sol(180): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(181): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(182): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(183): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(184): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(185): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(186): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(187): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(188): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(189): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(190): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(191): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,15) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(192): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(193): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(194): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(195): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(196): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(197): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(198): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(199): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(200): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(201): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(202): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(203): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(204): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(205): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(206): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(207): lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(6,2,11) <= 1
 sol(208): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(209): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(210): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(211): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(212): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(213): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(214): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(215): lambda(2,2,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(216): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(217): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(218): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(219): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(220): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(221): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(222): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(223): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(224): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(225): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(226): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(227): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(228): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(229): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,11) + lambda(2,3,13) <= 1
 sol(230): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(231): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,11) <= 1
 sol(232): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(233): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(2,3,13) <= 1
 sol(234): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(235): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(236): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(237): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(238): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(239): lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(240): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(241): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(242): lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(243): lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(244): lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(245): lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(246): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(247): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(248): lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(249): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(250): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(251): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(252): lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(253): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(254): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,4,14) + lambda(2,3,13) <= 1
 sol(255): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(256): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(257): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(258): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(259): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(260): lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(261): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(262): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(263): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(264): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(265): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(266): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(267): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(268): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(269): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(270): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(271): lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(272): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(273): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(274): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(275): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(276): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(277): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(278): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(279): lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(280): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(281): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(282): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(283): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(284): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(285): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13) + lambda(5,3,13)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(286): lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(287): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13) + lambda(5,3,13)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(288): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(289): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13) + lambda(5,3,13)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(290): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(291): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(292): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(293): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(294): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(295): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(296): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(297): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(298): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(299): lambda(6,1,15) + lambda(6,2,8) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(300): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(301): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(302): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(303): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(304): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(305): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(306): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(1,2,15) <= 1
 sol(307): lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(308): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(309): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(310): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(311): lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(312): lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(313): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(314): lambda(4,1,13) + lambda(6,2,8) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(315): lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(316): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(317): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(318): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,10) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(319): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(320): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(321): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(322): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(323): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(324): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(325): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(326): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,9) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(327): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(328): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   <= 1
 sol(329): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(330): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(331): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(332): lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(333): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(6,2,11) <= 1
 sol(334): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(335): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(336): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(5,2,13) <= 1
 sol(337): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(338): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   <= 1
 sol(339): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15)
   <= 1
 sol(340): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(341): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   <= 1
 sol(342): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(343): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(344): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(345): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(346): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(347): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(348): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(349): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(350): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13)
   <= 1
 sol(351): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(352): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(353): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(354): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(5,2,13) <= 1
 sol(355): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,2,13) <= 1
 sol(356): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(6,2,11) <= 1
 sol(357): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(358): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(1,2,15) <= 1
 sol(359): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(360): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(1,2,15) <= 1
 sol(361): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(5,2,13) <= 1
 sol(362): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(5,2,13) <= 1
 sol(363): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(6,2,11)
   <= 1
 sol(364): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(5,2,13) <= 1
 sol(365): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(366): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(367): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(6,2,11)
   <= 1
 sol(368): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(369): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(370): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(371): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13)
   <= 1
 sol(372): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13)
   <= 1
 sol(373): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13)
   <= 1
 sol(374): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(375): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   <= 1
 sol(376): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(377): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   <= 1
 sol(378): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(379): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(380): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(381): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(382): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(1,2,15) <= 1
 sol(383): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(5,2,13) <= 1
 sol(384): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   <= 1
 sol(385): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(6,2,11)
   <= 1
 sol(386): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) <= 1
 sol(387): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) <= 1
 sol(388): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(389): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(390): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(4,2,13) <= 1
 sol(391): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(392): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(393): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(394): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(395): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(396): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,10) <= 1
 sol(397): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(398): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(399): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(400): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(401): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(402): lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(403): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(404): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(405): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(406): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(407): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,3,10)
   <= 1
 sol(408): lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(409): lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(410): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,3,10)
   <= 1
 sol(411): lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(412): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(413): lambda(2,2,8) + lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(414): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(415): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(416): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(417): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(418): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(419): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(420): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(421): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(422): lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(423): lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(424): lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(425): lambda(4,1,13) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(426): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(427): lambda(4,1,13) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(428): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14) + lambda(1,3,10) <= 1
 sol(429): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(430): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(431): lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(432): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(433): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(434): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(435): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(436): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(437): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(438): lambda(3,3,11) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(439): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(440): lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(441): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(442): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(443): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(444): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(1,2,15) <= 1
 sol(445): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(446): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,13) <= 1
 sol(447): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(448): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(449): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(450): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14) + lambda(1,3,10) <= 1
 sol(451): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(452): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(453): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(454): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(455): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(456): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(457): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(458): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(459): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(460): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(461): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14) <= 1
 sol(462): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) <= 1
 sol(463): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) <= 1
 sol(464): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) <= 1
 sol(465): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(466): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) <= 1
 sol(467): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(468): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(469): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   <= 1
 sol(470): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(471): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   <= 1
 sol(472): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) <= 1
 sol(473): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(474): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(475): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(476): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(477): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(478): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(479): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(1,2,15) <= 1
 sol(480): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(481): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,2,13) <= 1
 sol(482): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,2,13) <= 1
 sol(483): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,2,13) <= 1
 sol(484): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(485): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(486): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(487): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(488): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(489): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(490): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(491): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(492): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(493): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) <= 1
 sol(494): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13) <= 1
 sol(495): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) <= 1
 sol(496): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13) <= 1
 sol(497): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(498): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(499): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(1,2,15) <= 1
 sol(500): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,2,13) <= 1
 sol(501): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,2,13) <= 1
 sol(502): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(503): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(504): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(505): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(506): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(507): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(6,6,14) + lambda(3,4,11) + lambda(6,2,11)
   <= 1
 sol(508): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(509): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(5,3,15) <= 1
 sol(510): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(511): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(512): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(513): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(514): lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(6,2,11) <= 1
 sol(515): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(516): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(517): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(518): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(519): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(520): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14) <= 1
 sol(521): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) <= 1
 sol(522): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) <= 1
 sol(523): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(524): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(525): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(526): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(527): lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(528): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(529): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(530): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(531): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(532): lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(533): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) <= 1
 sol(534): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(535): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(536): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,2,13) <= 1
 sol(537): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(538): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(539): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(540): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(541): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(542): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(543): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(544): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(545): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(546): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(547): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(548): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,2,13) <= 1
 sol(549): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(550): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(551): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(552): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   <= 1
 sol(553): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(554): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(555): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(556): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(557): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(558): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(559): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(560): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(561): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(562): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(563): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(564): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(565): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(566): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(567): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(568): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(569): lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(570): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(571): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(572): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(573): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(574): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(575): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(576): lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(577): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(578): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(579): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(580): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(581): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(582): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(583): lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(584): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(585): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(586): lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(587): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(588): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(589): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(590): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(591): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(592): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(593): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(594): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(595): lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(596): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(597): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(598): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(599): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(600): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(601): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(602): lambda(2,2,8) + lambda(1,1,8) + lambda(3,3,13) + lambda(6,6,9)
   + lambda(1,1,14) <= 1
 sol(603): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(604): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(605): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(606): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(607): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(608): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(609): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(610): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(611): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(612): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(613): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(614): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(615): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(616): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(617): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(618): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(619): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(620): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(621): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(622): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(623): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) <= 1
 sol(624): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) <= 1
 sol(625): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) <= 1
 sol(626): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15)
   <= 1
 sol(627): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) <= 1
 sol(628): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(629): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(630): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(631): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(632): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(633): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(634): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(635): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(636): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,15) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(637): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(638): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(639): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(640): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(641): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(642): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(643): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(644): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(645): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(646): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(647): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) <= 1
 sol(648): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(649): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(650): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(651): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(652): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(653): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(654): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(655): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(656): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(657): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(658): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(659): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(660): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(661): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(662): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(663): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(664): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(665): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(666): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(667): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(668): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(669): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(670): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,4,14) <= 1
 sol(671): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(672): lambda(4,1,13) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(673): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(674): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(675): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(676): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(677): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(678): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(679): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(680): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(681): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(682): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(683): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(684): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(685): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(686): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(687): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(688): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(689): lambda(4,1,13) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(690): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(691): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(692): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(693): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(694): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(695): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) <= 1
 sol(696): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(697): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(698): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(699): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(700): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(701): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(702): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(703): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(704): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(705): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(706): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(707): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(708): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(709): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,2,13)
   <= 1
 sol(710): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,2,13)
   <= 1
 sol(711): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,2,13)
   <= 1
 sol(712): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(1,2,15)
   <= 1
 sol(713): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(714): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,2,13)
   <= 1
 sol(715): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,2,13)
   <= 1
 sol(716): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(1,2,15)
   <= 1
 sol(717): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(718): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(719): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(720): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(721): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,2,13)
   <= 1
 sol(722): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(723): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(724): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(725): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(726): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(727): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11)
   <= 1
 sol(728): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(729): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(6,6,14) <= 1
 sol(730): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(731): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(732): lambda(6,1,15) + lambda(3,3,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(733): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(734): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(735): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(736): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(737): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(738): lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(739): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(740): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(741): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(742): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15)
   <= 1
 sol(743): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(744): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(745): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(746): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(747): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(5,2,13) <= 1
 sol(748): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(749): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(750): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(751): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(752): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(753): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(754): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(755): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(756): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(757): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(758): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(759): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(760): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(761): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(762): lambda(1,1,8) + lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(763): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(764): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(765): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(766): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(767): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(768): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(769): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(770): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(771): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(6,2,11) <= 1
 sol(772): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(773): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(774): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(775): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(776): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(777): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(778): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(779): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) <= 1
 sol(780): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(781): lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(782): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(783): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(784): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(785): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(786): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(787): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(788): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(789): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(790): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(791): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) <= 1
 sol(792): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(793): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(794): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(795): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(796): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(797): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(798): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(799): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(800): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(801): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(802): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(803): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(804): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(805): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(806): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(807): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(808): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(809): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(810): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(811): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(812): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(813): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(814): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(815): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(816): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(817): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(818): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(819): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(820): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(821): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(822): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(823): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(824): lambda(1,1,8) + lambda(6,1,15) + lambda(3,3,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(825): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(826): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(827): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(828): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(829): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(830): lambda(1,1,8) + lambda(6,1,15) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(831): lambda(1,1,8) + lambda(1,1,14) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(832): lambda(1,1,8) + lambda(1,1,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(833): lambda(1,1,8) + lambda(1,1,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(834): lambda(1,1,8) + lambda(4,1,13) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(835): lambda(1,1,8) + lambda(4,1,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(836): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(837): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(838): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(839): lambda(1,1,8) + lambda(4,1,13) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(840): lambda(1,1,8) + lambda(4,1,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(841): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(842): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(843): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(844): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(845): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(846): lambda(1,1,8) + lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(847): lambda(1,1,8) + lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(848): lambda(1,1,8) + lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(849): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(850): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(851): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(852): lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(853): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(854): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(855): lambda(1,1,8) + lambda(6,6,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(856): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(857): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(858): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(859): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(860): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(1,2,15)
   <= 1
 sol(861): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(862): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(863): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(864): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(865): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(866): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(867): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(868): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(869): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(870): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(871): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(872): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(873): lambda(1,1,8) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(874): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(875): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(876): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(877): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(878): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(879): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(880): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(881): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(882): lambda(1,1,8) + lambda(6,6,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(883): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(884): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(885): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(886): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(887): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(888): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(889): lambda(6,4,11) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(890): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(891): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(892): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(893): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(894): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(895): lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(896): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(897): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(898): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(899): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(900): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(901): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(902): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(903): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(904): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(905): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(906): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(907): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(908): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(909): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(910): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(911): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15) + lambda(1,3,10) <= 1
 sol(912): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(913): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(914): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(915): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(916): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(917): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(918): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(919): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(920): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(921): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(922): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(923): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(924): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(925): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(926): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(927): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(928): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,10) <= 1
 sol(929): lambda(6,4,11) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(930): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(931): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(932): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(933): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(934): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(935): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(936): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(937): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(938): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(939): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(940): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(941): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(942): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(943): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(944): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(945): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(946): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(947): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(948): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(949): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(4,2,13) <= 1
 sol(950): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(951): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(952): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(953): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(954): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(955): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(956): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(957): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(958): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(959): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(960): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(961): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(1,2,15) <= 1
 sol(962): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(963): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(964): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(965): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,6,14) <= 1
 sol(966): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(967): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(968): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(969): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(970): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(971): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(972): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(973): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(974): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(975): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(976): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(977): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(978): lambda(1,1,8) + lambda(1,1,14) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(979): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(980): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(981): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(982): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(983): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(984): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(985): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(986): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(987): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(988): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(989): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(990): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(991): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(992): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(993): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(994): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13) + lambda(1,3,10) <= 1
 sol(995): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(996): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(997): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(998): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(999): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(1000): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1001): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(2,3,8) <= 1
 sol(1002): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1003): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1004): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1005): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1006): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1007): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) <= 1
 sol(1008): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1009): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,2,13) <= 1
 sol(1010): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1011): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1012): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1013): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1014): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(1015): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1016): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) <= 1
 sol(1017): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1018): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1019): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1020): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1021): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(1022): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1023): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1024): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1025): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14) <= 1
 sol(1026): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1027): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1028): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1029): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1030): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1031): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1032): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1033): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) <= 1
 sol(1034): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1035): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(1036): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1037): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1038): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1039): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1040): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1041): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1042): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(1043): lambda(1,1,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1044): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1045): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1046): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1047): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1048): lambda(2,2,8) + lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1049): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1050): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1051): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(1052): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1053): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1054): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1055): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(1056): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1057): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(1058): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1059): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1060): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(1061): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(1062): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1063): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1064): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1065): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1066): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1067): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1068): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(1,2,15) <= 1
 sol(1069): lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1070): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(1071): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1072): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1073): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1074): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1075): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1076): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1077): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1078): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(1,3,8) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1079): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,3,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1080): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1081): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,3,10) <= 1
 sol(1082): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1083): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1084): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1085): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1086): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(1087): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1088): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1089): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1090): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1091): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1092): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1093): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1094): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1095): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1096): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(1097): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1098): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1099): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1100): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1101): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1102): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1103): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1104): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1105): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1106): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1107): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1108): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(1,3,10) <= 1
 sol(1109): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1110): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(1111): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1112): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,3,13) <= 1
 sol(1113): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) <= 1
 sol(1114): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(1115): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1116): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1117): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1118): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1119): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1120): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(1121): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1122): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1123): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1124): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(1125): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1126): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1127): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(1128): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1129): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1130): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1131): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1132): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1133): lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(1134): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1135): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1136): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(6,2,11) <= 1
 sol(1137): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1138): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) <= 1
 sol(1139): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1140): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1141): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1142): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(1143): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1144): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1145): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1146): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1147): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1148): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1149): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1150): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1151): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(1152): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1153): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(1154): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1155): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1156): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1157): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1158): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1159): lambda(2,2,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(1160): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1161): lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15) <= 1
 sol(1162): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1163): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) <= 1
 sol(1164): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(1165): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(1166): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1167): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(1168): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1169): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(1170): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1171): lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(1172): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(1173): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1174): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1175): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1176): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1177): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1178): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1179): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1180): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1181): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1182): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(1183): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1184): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1185): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1186): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1187): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1188): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1189): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1190): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1191): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1192): lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1193): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1194): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1195): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) <= 1
 sol(1196): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1197): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1198): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1199): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1200): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1201): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1202): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1203): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1204): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1205): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1206): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1207): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1208): lambda(2,2,8) + lambda(1,1,8) + lambda(6,1,15) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(1209): lambda(6,1,15) + lambda(3,3,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1210): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1211): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1212): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,2,13) <= 1
 sol(1213): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1214): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1215): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1216): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1217): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1218): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(1219): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1220): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1221): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1222): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(1223): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1224): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1225): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1226): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1227): lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(1228): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1229): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1230): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(1231): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1232): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1233): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1234): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1235): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1236): lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1237): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(1238): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1239): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1240): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1241): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1242): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1243): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1244): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1245): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1246): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(1247): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) <= 1
 sol(1248): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1249): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1250): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1251): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1252): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1253): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1254): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1255): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1256): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1257): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(1258): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1259): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(6,6,14) <= 1
 sol(1260): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(1261): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1262): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1263): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1264): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(1265): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1266): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1267): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1268): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1269): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1270): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1271): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1272): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1273): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1274): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1275): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1276): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(1277): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1278): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(1279): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1280): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1281): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1282): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(1283): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1284): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1285): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1286): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1287): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1288): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1289): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1290): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1291): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(1292): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(1293): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(1294): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1295): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1296): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(1,2,15) <= 1
 sol(1297): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1298): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1299): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1300): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1301): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1302): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1303): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1304): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1305): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(1306): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1307): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1308): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(1309): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1310): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1311): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(1312): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1313): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1314): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) <= 1
 sol(1315): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1316): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1317): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1318): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1319): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1320): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(1321): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(1322): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1323): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(1324): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(1325): lambda(6,1,15) + lambda(3,3,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1326): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1327): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1328): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1329): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1330): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1331): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1332): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(1333): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1334): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(1335): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1336): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1337): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1338): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(1339): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1340): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1341): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1342): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1343): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1344): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1345): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1346): lambda(3,3,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1347): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1348): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1349): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,9) + lambda(2,3,13) <= 1
 sol(1350): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1351): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1352): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1353): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(4,2,13) <= 1
 sol(1354): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1355): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1356): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1357): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1358): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1359): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1360): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1361): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1362): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1363): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1364): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1365): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1366): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1367): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1368): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1369): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1370): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1371): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1372): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(1373): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1374): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1375): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1376): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1377): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1378): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(1379): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1380): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1381): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(1382): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1383): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1384): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1385): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1386): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1387): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1388): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1389): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1390): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1391): lambda(6,6,9) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1392): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1393): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1394): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1395): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1396): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1397): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(1398): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1399): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1400): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(1401): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1402): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1403): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1404): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1405): lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(1406): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1407): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1408): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(1409): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1410): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1411): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1412): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1413): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1414): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1415): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(1416): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1417): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1418): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1419): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1420): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1421): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1422): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1423): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1424): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1425): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1426): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1427): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1428): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1429): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1430): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1431): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1432): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1433): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1434): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1435): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1436): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1437): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1438): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1439): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1440): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1441): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1442): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1443): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1444): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(1445): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1446): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1447): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1448): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(1449): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1450): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1451): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(1452): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1453): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1454): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1455): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(6,2,11) <= 1
 sol(1456): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1457): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1458): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(1459): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1460): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1461): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(1462): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1463): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1464): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1465): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1466): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(1467): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1468): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1469): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1470): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1471): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(1472): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1473): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1474): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1475): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1476): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1477): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1478): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1479): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1480): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1481): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1482): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1483): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1484): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1485): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1486): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1487): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1488): lambda(1,1,14) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(1489): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(1490): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(4,2,13) <= 1
 sol(1491): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(1492): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1493): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1494): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(1495): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1496): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1497): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1498): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1499): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1500): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(5,3,13) <= 1
 sol(1501): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1502): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1503): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1504): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1505): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1506): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1507): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1508): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1509): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1510): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(1511): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) <= 1
 sol(1512): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1513): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(1514): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1515): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1516): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(1517): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1518): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1519): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1520): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1521): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1522): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1523): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1524): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1525): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1526): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1527): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1528): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) <= 1
 sol(1529): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1530): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(1531): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1532): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(1533): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(1,3,10) <= 1
 sol(1534): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1535): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1536): lambda(4,1,13) + lambda(6,6,9) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1537): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1538): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,10) <= 1
 sol(1539): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1540): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(1541): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1542): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(1543): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1544): lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1545): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1546): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(1,2,15) <= 1
 sol(1547): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1548): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1549): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1550): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1551): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1552): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1553): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1554): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1555): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1556): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1557): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1558): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1559): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1560): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) <= 1
 sol(1561): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1562): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1563): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1564): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(4,2,13) <= 1
 sol(1565): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1566): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15)
   <= 1
 sol(1567): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1568): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1569): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1570): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1571): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(1572): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1573): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,4,14) <= 1
 sol(1574): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(1,3,10) <= 1
 sol(1575): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1576): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1577): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1578): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1579): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1580): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1581): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1582): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1583): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1584): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1585): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1586): lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1587): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1588): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1589): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1590): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1591): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1592): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1593): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(1,2,15) <= 1
 sol(1594): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1595): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1596): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1597): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1598): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1599): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1600): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(1601): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1602): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1603): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(1604): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1605): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1606): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1607): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1608): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1609): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1610): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(1611): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1612): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1613): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(1614): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1615): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1616): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1617): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1618): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1619): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1620): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,9) + lambda(2,3,13) <= 1
 sol(1621): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1622): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1623): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1624): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1625): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1626): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1627): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1628): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1629): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1630): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1631): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1632): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1633): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(1634): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(1635): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1636): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) <= 1
 sol(1637): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1638): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1639): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1640): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1641): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1642): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(1643): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1644): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1645): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1646): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1647): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14)
   <= 1
 sol(1648): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1649): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1650): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1651): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1652): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1653): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(6,2,11) <= 1
 sol(1654): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1655): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1656): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(1657): lambda(1,1,8) + lambda(6,1,15) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1658): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1659): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1660): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(1661): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1662): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1663): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1664): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(1665): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1666): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1667): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1668): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(1669): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1670): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1671): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1672): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1673): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1674): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1675): lambda(2,2,8) + lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1676): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1677): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1678): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1679): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(1680): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1681): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1682): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(1683): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1684): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1685): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(1686): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1687): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1688): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1689): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1690): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1691): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1692): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(1693): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1694): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(1695): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(1696): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1697): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1698): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1699): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1700): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1701): lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1702): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1703): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1704): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1705): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1706): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1707): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1708): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1709): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1710): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(1711): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(1712): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1713): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1714): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) <= 1
 sol(1715): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1716): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(1717): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1718): lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(2,3,8) <= 1
 sol(1719): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1720): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1721): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1722): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1723): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1724): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1725): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(1,3,8) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1726): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1727): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1728): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1729): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1730): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1731): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1732): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1733): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(1,2,15) <= 1
 sol(1734): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1735): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1736): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1737): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(1738): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1739): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1740): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1741): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1742): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1743): lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1744): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1745): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1746): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(1747): lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(1748): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1749): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1750): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1751): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1752): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1753): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1754): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1755): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1756): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1757): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(1758): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1759): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1760): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1761): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(1762): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1763): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1764): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(1765): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14)
   <= 1
 sol(1766): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1767): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(1768): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1769): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1770): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1771): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1772): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1773): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1774): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1775): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(1776): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(6,2,11)
   <= 1
 sol(1777): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1778): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(1779): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(6,2,11) <= 1
 sol(1780): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1781): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(1782): lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(1783): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1784): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1785): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1786): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1787): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1788): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(1789): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(1790): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1791): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1792): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1793): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(1,2,15) <= 1
 sol(1794): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(1795): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(1796): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1797): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(1798): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(5,2,13) <= 1
 sol(1799): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1800): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(1801): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) <= 1
 sol(1802): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1803): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1804): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1805): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1806): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1807): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1808): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1809): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) <= 1
 sol(1810): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1811): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(1812): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1813): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(1814): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1815): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1816): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1817): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1818): lambda(1,1,8) + lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(2,3,8) <= 1
 sol(1819): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13)
   <= 1
 sol(1820): lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(1821): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1822): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1823): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1824): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1825): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1826): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1827): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) <= 1
 sol(1828): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(1829): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1830): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1831): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(1832): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1833): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1834): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1835): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1836): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1837): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1838): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1839): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1840): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(1841): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1842): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1843): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1844): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1845): lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1846): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1847): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1848): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1849): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1850): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1851): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1852): lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1853): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(1854): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1855): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(1856): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1857): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(1858): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1859): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1860): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1861): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1862): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(4,3,14) + lambda(3,3,9) <= 1
 sol(1863): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(1864): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1865): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1866): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1867): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1868): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1869): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(1870): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1871): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1872): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1873): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1874): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1875): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1876): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1877): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1878): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(1879): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1880): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1881): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1882): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1883): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(1884): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(1885): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1886): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1887): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1888): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1889): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1890): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1891): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1892): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1893): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(1894): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1895): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(1896): lambda(1,1,8) + lambda(6,1,15) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(1897): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(1898): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1899): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(1900): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1901): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(1902): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1903): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1904): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1905): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1906): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1907): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(1908): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1909): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1910): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1911): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(1912): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1913): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(2,3,8) <= 1
 sol(1914): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1915): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1916): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1917): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(1918): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1919): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1920): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1921): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(1922): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(1923): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(1924): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1925): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(1926): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1927): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1928): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1929): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(1930): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1931): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1932): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(1933): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(1934): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1935): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1936): lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(1937): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1938): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1939): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(1940): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1941): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1942): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(1943): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1944): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(1945): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1946): lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(1947): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1948): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(1949): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1950): lambda(6,4,11) + lambda(1,3,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(1951): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1952): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(1953): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1954): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1955): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1956): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(1957): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1958): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1959): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,4,14) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(1960): lambda(3,3,13) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(1961): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(1962): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(1963): lambda(1,1,8) + lambda(6,1,15) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(1964): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1965): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(1966): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(1967): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(1968): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(1969): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1970): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1971): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1972): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1973): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(6,2,11) <= 1
 sol(1974): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1975): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(1976): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(1977): lambda(4,1,13) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(1978): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(1979): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(1980): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(1981): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1982): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1983): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1984): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(1985): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(1986): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(2,3,13) <= 1
 sol(1987): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1988): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(2,3,8) <= 1
 sol(1989): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(1990): lambda(1,1,8) + lambda(6,1,15) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) <= 1
 sol(1991): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(1992): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(1993): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(1994): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(1995): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(1996): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(1997): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(1998): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,2,13) <= 1
 sol(1999): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2000): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2001): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2002): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(5,3,13) <= 1
 sol(2003): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2004): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2005): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2006): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2007): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2008): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2009): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2010): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) <= 1
 sol(2011): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2012): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2013): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2014): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(2015): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2016): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2017): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2018): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2019): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2020): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2021): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2022): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2023): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2024): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2025): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2026): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2027): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2028): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2029): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(2030): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2031): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2032): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2033): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2034): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(2035): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2036): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2037): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2038): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2039): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2040): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(2041): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2042): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2043): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2044): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2045): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2046): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(2047): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13)
   <= 1
 sol(2048): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2049): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2050): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2051): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,11) + lambda(5,3,15) <= 1
 sol(2052): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2053): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2054): lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2055): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2056): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(1,2,15) <= 1
 sol(2057): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2058): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2059): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2060): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(2061): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(2062): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2063): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2064): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2065): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2066): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2067): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2068): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(2069): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2070): lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(2071): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2072): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2073): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2074): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2075): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2076): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(6,2,11) <= 1
 sol(2077): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2078): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2079): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2080): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2081): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2082): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2083): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2084): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2085): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2086): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2087): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(3,4,11) <= 1
 sol(2088): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(6,2,11) <= 1
 sol(2089): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2090): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,2,13) <= 1
 sol(2091): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2092): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2093): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2094): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2095): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2096): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2097): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2098): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2099): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,3,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2100): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(1,2,15) <= 1
 sol(2101): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2102): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2103): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2104): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2105): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2106): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2107): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2108): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2109): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2110): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2111): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(2,3,13) <= 1
 sol(2112): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2113): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2114): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2115): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2116): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2117): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2118): lambda(1,1,8) + lambda(6,1,15) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2119): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) <= 1
 sol(2120): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2121): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2122): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2123): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) <= 1
 sol(2124): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2125): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2126): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2127): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2128): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2129): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(2130): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2131): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2132): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2133): lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(2134): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2135): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2136): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2137): lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15)
   <= 1
 sol(2138): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2139): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2140): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2141): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(2142): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2143): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2144): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2145): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2146): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2147): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2148): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(2149): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2150): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2151): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2152): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2153): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2154): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2155): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(2156): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10)
   <= 1
 sol(2157): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2158): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2159): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2160): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2161): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2162): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2163): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2164): lambda(6,4,11) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2165): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(2166): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2167): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2168): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2169): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2170): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(2,3,8) <= 1
 sol(2171): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2172): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2173): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2174): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2175): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2176): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2177): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2178): lambda(6,4,11) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2179): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(2180): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,2,13) <= 1
 sol(2181): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2182): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2183): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,4,11) <= 1
 sol(2184): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2185): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2186): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2187): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2188): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2189): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2190): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2191): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2192): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2193): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(2194): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(2195): lambda(3,3,13) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(2196): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2197): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2198): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,3,10) <= 1
 sol(2199): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2200): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2201): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2202): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2203): lambda(4,1,13) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2204): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2205): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2206): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2207): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2208): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2209): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(5,3,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2210): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2211): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2212): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2213): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2214): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2215): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2216): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(2217): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2218): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2219): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2220): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(2221): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2222): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2223): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2224): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2225): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2226): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2227): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2228): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2229): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2230): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(2231): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(2232): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(2233): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2234): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2235): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2236): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2237): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2238): lambda(6,4,11) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2239): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2240): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2241): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(2242): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2243): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,11) <= 1
 sol(2244): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2245): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,6,9)
   + lambda(1,1,14) <= 1
 sol(2246): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2247): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2248): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2249): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2250): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2251): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2252): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2253): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(2254): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2255): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(2256): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(2257): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2258): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2259): lambda(3,3,11) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2260): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2261): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2262): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,13) <= 1
 sol(2263): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2264): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(2265): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2266): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2267): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2268): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2269): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2270): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2271): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2272): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2273): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2274): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2275): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2276): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2277): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(2278): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2279): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2280): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2281): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(2282): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2283): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2284): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2285): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2286): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2287): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2288): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2289): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2290): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2291): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2292): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2293): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2294): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2295): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2296): lambda(2,2,8) + lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) <= 1
 sol(2297): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,3,10) <= 1
 sol(2298): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(2299): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2300): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2301): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2302): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2303): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2304): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2305): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2306): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2307): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) <= 1
 sol(2308): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(2309): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2310): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2311): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(2312): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(2313): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2314): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2315): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(2316): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2317): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2318): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(2319): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2320): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2321): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2322): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2323): lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2324): lambda(1,1,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2325): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(2326): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2327): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2328): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2329): lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(2330): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2331): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2332): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2333): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2334): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(6,2,11) <= 1
 sol(2335): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2336): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(2337): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2338): lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2339): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2340): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(6,2,11) <= 1
 sol(2341): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2342): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2343): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2344): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2345): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(2346): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2347): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2348): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2349): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(2350): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2351): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2352): lambda(6,4,11) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2353): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(2354): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(2355): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2356): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2357): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2358): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2359): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(2360): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(2361): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(2362): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2363): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2364): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2365): lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(2366): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13)
   <= 1
 sol(2367): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2368): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(2369): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(2370): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2371): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2372): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2373): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2374): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2375): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2376): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(2377): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) <= 1
 sol(2378): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2379): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2380): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2381): lambda(2,2,8) + lambda(1,1,8) + lambda(6,1,15) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(2382): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2383): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2384): lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(2385): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2386): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2387): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2388): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2389): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2390): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2391): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2392): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2393): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2394): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(2395): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) <= 1
 sol(2396): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2397): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,6,10) + lambda(2,3,13) <= 1
 sol(2398): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2399): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2400): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2401): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2402): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2403): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2404): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2405): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2406): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2407): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2408): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2409): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(2410): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(2411): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(2412): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(2413): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(2414): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2415): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2416): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2417): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2418): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2419): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(2420): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2421): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(2422): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2423): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13) + lambda(3,6,9)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8) <= 1
 sol(2424): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2425): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2426): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2427): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2428): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2429): lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2430): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2431): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2432): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(1,3,10) <= 1
 sol(2433): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2434): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2435): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2436): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(6,2,11) <= 1
 sol(2437): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2438): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2439): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) <= 1
 sol(2440): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2441): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2442): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2443): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2444): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2445): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2446): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2447): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2448): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2449): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2450): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2451): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2452): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2453): lambda(2,2,8) + lambda(4,1,13) + lambda(1,3,8) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2454): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2455): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2456): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2457): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2458): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2459): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2460): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2461): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2462): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2463): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2464): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2465): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2466): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2467): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2468): lambda(2,2,8) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(2469): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2470): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2471): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) <= 1
 sol(2472): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2473): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2474): lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2475): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2476): lambda(1,1,8) + lambda(1,1,14) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2477): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2478): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2479): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2480): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2481): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2482): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2483): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2484): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2485): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(2486): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2487): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2488): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2489): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2490): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2491): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(2492): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2493): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2494): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2495): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2496): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2497): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2498): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2499): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2500): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2501): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2502): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2503): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2504): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2505): lambda(3,3,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2506): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2507): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2508): lambda(3,3,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2509): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2510): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(2511): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2512): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2513): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(2514): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2515): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2516): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2517): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,11) + lambda(5,3,13) <= 1
 sol(2518): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2519): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2520): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(2521): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2522): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2523): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2524): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2525): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(2526): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2527): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2528): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2529): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2530): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2531): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2532): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2533): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2534): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2535): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2536): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(2537): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) <= 1
 sol(2538): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(6,2,11) <= 1
 sol(2539): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2540): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2541): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2542): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2543): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2544): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2545): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2546): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2547): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2548): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(2549): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2550): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2551): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,13) <= 1
 sol(2552): lambda(3,3,13) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2553): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2554): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2555): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2556): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2557): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2558): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2559): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2560): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2561): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2562): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2563): lambda(2,2,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,4,11) + lambda(4,3,14) <= 1
 sol(2564): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2565): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2566): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) <= 1
 sol(2567): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(2568): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2569): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2570): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2571): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(2572): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2573): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2574): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2575): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2576): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2577): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(2578): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(2579): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2580): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2581): lambda(6,4,11) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2582): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2583): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2584): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2585): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2586): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9)
   + lambda(3,4,14) + lambda(5,3,13) <= 1
 sol(2587): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2588): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2589): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2590): lambda(6,1,15) + lambda(1,3,8) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) + lambda(1,3,10)
   <= 1
 sol(2591): lambda(6,4,11) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2592): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(3,3,9) <= 1
 sol(2593): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(4,1,8) <= 1
 sol(2594): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2595): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2596): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2597): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2598): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2599): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2600): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2601): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2602): lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2603): lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2604): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2605): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2606): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2607): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(2608): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2609): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2610): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2611): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2612): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2613): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2614): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2615): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2616): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2617): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2618): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2619): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2620): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2621): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2622): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2623): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2624): lambda(3,3,13) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2625): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2626): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2627): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2628): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(3,4,11)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2629): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2630): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2631): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2632): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2633): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2634): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2635): lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2636): lambda(1,1,8) + lambda(4,1,13) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2637): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2638): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(1,2,15) <= 1
 sol(2639): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2640): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(4,2,13) <= 1
 sol(2641): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2642): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2643): lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2644): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2645): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2646): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2647): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(6,2,11) <= 1
 sol(2648): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2649): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2650): lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(2651): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2652): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2653): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2654): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2655): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2656): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2657): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2658): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2659): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2660): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2661): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2662): lambda(2,2,8) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,14)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(2663): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2664): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2665): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2666): lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(2667): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2668): lambda(1,1,8) + lambda(6,6,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(4,2,8) + lambda(4,3,14) + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2669): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2670): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2671): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2672): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2673): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2674): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2675): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(2676): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2677): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(2678): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2679): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2680): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(2681): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,2,13) <= 1
 sol(2682): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(2683): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2684): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) <= 1
 sol(2685): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2686): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2687): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(3,4,14) + lambda(2,3,13) <= 1
 sol(2688): lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15) <= 1
 sol(2689): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(2690): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,11) <= 1
 sol(2691): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2692): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2693): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(2694): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,13) <= 1
 sol(2695): lambda(2,2,8) + lambda(3,3,11) + lambda(6,4,11)
   + lambda(4,1,13) + lambda(1,3,8) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2696): lambda(2,2,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(3,3,9) <= 1
 sol(2697): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(2698): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2699): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2700): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2701): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2702): lambda(3,3,13) + lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(2703): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2704): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2705): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(1,2,15) <= 1
 sol(2706): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2707): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2708): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2709): lambda(2,2,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,6,14) <= 1
 sol(2710): lambda(6,1,15) + lambda(1,1,9) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2711): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2712): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2713): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2714): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2715): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2716): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2717): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2718): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2719): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2720): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2721): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) <= 1
 sol(2722): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2723): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2724): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,3,14) <= 1
 sol(2725): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) <= 1
 sol(2726): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2727): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2728): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2729): lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2730): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2731): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2732): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(2733): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2734): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2735): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) <= 1
 sol(2736): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(2737): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2738): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(2739): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(2740): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2741): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2742): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2743): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2744): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2745): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2746): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2747): lambda(1,1,8) + lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(2748): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2749): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2750): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2751): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14)
   + lambda(6,2,11) <= 1
 sol(2752): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2753): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2754): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2755): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8) <= 1
 sol(2756): lambda(6,4,11) + lambda(1,1,9) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2757): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(2758): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(2759): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2760): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(1,2,15) <= 1
 sol(2761): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2762): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2763): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2764): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(6,2,11)
   <= 1
 sol(2765): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2766): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2767): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(2768): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(3,4,14)
   + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2769): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2770): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2771): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(2772): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2773): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(2774): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2775): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) <= 1
 sol(2776): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2777): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2778): lambda(4,1,13) + lambda(6,6,8) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2779): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2780): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2781): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2782): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2783): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2784): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,10) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(2785): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2786): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(2787): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2788): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2789): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2790): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2791): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(2792): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2793): lambda(1,1,8) + lambda(1,1,14) + lambda(6,6,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(2794): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2795): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11)
   <= 1
 sol(2796): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2797): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2798): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(2799): lambda(3,3,11) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(2800): lambda(6,4,11) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2801): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2802): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2803): lambda(1,1,14) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(2804): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(2805): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2806): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2807): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2808): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(2809): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(2810): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2811): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2812): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(1,2,15) <= 1
 sol(2813): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2814): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(2815): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(2816): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2817): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2818): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8) <= 1
 sol(2819): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2820): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2821): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2822): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2823): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(2824): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2825): lambda(4,1,13) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2826): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2827): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(1,3,10) <= 1
 sol(2828): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2829): lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2830): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2831): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2832): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2833): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2834): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(4,2,13) <= 1
 sol(2835): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(2836): lambda(1,1,8) + lambda(4,1,13) + lambda(2,3,11) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(5,2,13) <= 1
 sol(2837): lambda(6,4,11) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2838): lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2839): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(2840): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(2841): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(6,6,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2842): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2843): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2844): lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2845): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2846): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2847): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2848): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2849): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(2850): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2851): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2852): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2853): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(2854): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2855): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2856): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2857): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2858): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) <= 1
 sol(2859): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2860): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2861): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2862): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2863): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(2,3,8) <= 1
 sol(2864): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(2865): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2866): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2867): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,13)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2868): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(2869): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(2870): lambda(3,3,13) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(4,2,8) + lambda(2,3,8) <= 1
 sol(2871): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2872): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(2,3,13) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2873): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2874): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2875): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2876): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(2877): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,2,13)
   <= 1
 sol(2878): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2879): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2880): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(2881): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(2882): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(2883): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2884): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2885): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8)
   + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(2886): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2887): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(2888): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(3,3,9) <= 1
 sol(2889): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,3,9) <= 1
 sol(2890): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(2891): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2892): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(2893): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) <= 1
 sol(2894): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(1,2,15) <= 1
 sol(2895): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2896): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2897): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2898): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2899): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2900): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2901): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2902): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2903): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2904): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2905): lambda(3,3,11) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2906): lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2907): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(2908): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(2909): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2910): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,3,10) <= 1
 sol(2911): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(2912): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15) <= 1
 sol(2913): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2914): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(2915): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2916): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2917): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2918): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(2919): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2920): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(6,2,11) <= 1
 sol(2921): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(2922): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2923): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2924): lambda(3,3,11) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(2925): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(2926): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2927): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(2928): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(2929): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(6,4,14) <= 1
 sol(2930): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2931): lambda(6,1,15) + lambda(1,3,8) + lambda(4,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,2,13) + lambda(1,3,10)
   <= 1
 sol(2932): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,2,13) <= 1
 sol(2933): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2934): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(2935): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(2936): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(6,6,14) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2937): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2938): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(2939): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(5,2,13) <= 1
 sol(2940): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2941): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,4,14) + lambda(6,6,14)
   + lambda(5,2,13) <= 1
 sol(2942): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,4,14) + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(2943): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2944): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2945): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2946): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(2947): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(2948): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(2949): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(4,2,13) <= 1
 sol(2950): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(4,2,13) <= 1
 sol(2951): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(2952): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2953): lambda(3,3,11) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(2954): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(2955): lambda(3,3,11) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(1,2,15) <= 1
 sol(2956): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(5,2,13)
   <= 1
 sol(2957): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2958): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(2959): lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(2960): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2961): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2962): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2963): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(6,2,11) <= 1
 sol(2964): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(2965): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(6,2,11) <= 1
 sol(2966): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(2967): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(6,2,11) <= 1
 sol(2968): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2969): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) <= 1
 sol(2970): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) <= 1
 sol(2971): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) + lambda(3,3,9) <= 1
 sol(2972): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2973): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(2974): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2975): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2976): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,2,13) <= 1
 sol(2977): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,3,14) <= 1
 sol(2978): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(6,2,13) + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) <= 1
 sol(2979): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2980): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,9) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2981): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(2982): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(2983): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(6,2,11) <= 1
 sol(2984): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(6,2,11) <= 1
 sol(2985): lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2986): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(2987): lambda(1,1,14) + lambda(6,6,8) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(2988): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2989): lambda(2,2,8) + lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(1,1,14) + lambda(6,4,14) <= 1
 sol(2990): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(2991): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(2992): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(2993): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(2994): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14)
   + lambda(5,2,13) <= 1
 sol(2995): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(5,3,15) + lambda(4,2,13) <= 1
 sol(2996): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(4,2,13) <= 1
 sol(2997): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(5,2,13) <= 1
 sol(2998): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(2999): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(3000): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(3001): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(3002): lambda(6,4,11) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(3003): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(3004): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(3005): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(3006): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3007): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(4,3,14) + lambda(1,2,15) + lambda(6,2,11)
   + lambda(3,3,9) <= 1
 sol(3008): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3009): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(3010): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(3011): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3012): lambda(3,3,11) + lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,6,14) + lambda(6,2,11) <= 1
 sol(3013): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(3014): lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(3015): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(6,2,11) <= 1
 sol(3016): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(3017): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(3018): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3019): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(6,6,8)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(3020): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3021): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(3022): lambda(2,2,8) + lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(3023): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3024): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(3025): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3026): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,6,11) + lambda(4,3,14) <= 1
 sol(3027): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(3028): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,11) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(3029): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(3030): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,2,13) <= 1
 sol(3031): lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(2,3,13) + lambda(1,2,15) + lambda(6,2,11)
   <= 1
 sol(3032): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(3033): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3034): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   + lambda(4,2,13) <= 1
 sol(3035): lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) <= 1
 sol(3036): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(3037): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(4,2,13) <= 1
 sol(3038): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(3039): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(1,2,15) <= 1
 sol(3040): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(5,2,13) <= 1
 sol(3041): lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) <= 1
 sol(3042): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3043): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3044): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3045): lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3046): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3047): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3048): lambda(3,3,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3049): lambda(6,1,15) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3050): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3051): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3052): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,3,11) + lambda(3,6,9) + lambda(5,2,13) <= 1
 sol(3053): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3054): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3055): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3056): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3057): lambda(3,3,11) + lambda(6,4,11) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(1,2,15) <= 1
 sol(3058): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3059): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3060): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3061): lambda(1,1,8) + lambda(4,1,13) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(3062): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3063): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3064): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3065): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3066): lambda(2,2,8) + lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) <= 1
 sol(3067): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,4,14) + lambda(5,2,13) <= 1
 sol(3068): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3069): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3070): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3071): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(1,2,15) <= 1
 sol(3072): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(3073): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(5,3,13) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3074): lambda(4,1,13) + lambda(6,6,9) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(3075): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3076): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3077): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(3078): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3079): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(1,2,15) <= 1
 sol(3080): lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(6,2,11) <= 1
 sol(3081): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,9) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(1,2,15) <= 1
 sol(3082): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3083): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3084): lambda(1,1,8) + lambda(1,1,14) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3085): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3086): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3087): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3088): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(1,2,15) <= 1
 sol(3089): lambda(6,4,11) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(3090): lambda(3,3,13) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3091): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(3,6,11) + lambda(2,3,8) + lambda(2,3,13) <= 1
 sol(3092): lambda(6,2,8) + lambda(2,5,13) + lambda(3,4,14)
   + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3093): lambda(6,4,11) + lambda(6,1,15) + lambda(3,3,13)
   + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,2,13) <= 1
 sol(3094): lambda(3,3,11) + lambda(1,1,8) + lambda(6,4,11)
   + lambda(1,1,14) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(4,2,13) <= 1
 sol(3095): lambda(6,2,8) + lambda(2,5,13) + lambda(3,6,15)
   + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3096): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3097): lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3098): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3099): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,13) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(2,3,13) + lambda(1,3,10) <= 1
 sol(3100): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(3,3,9) <= 1
 sol(3101): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3102): lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3103): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11)
   + lambda(1,2,15) <= 1
 sol(3104): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3105): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3106): lambda(6,6,9) + lambda(1,1,14) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(4,2,13) <= 1
 sol(3107): lambda(2,2,8) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(5,3,13) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(3108): lambda(6,1,15) + lambda(1,3,8) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(1,2,15)
   + lambda(1,3,10) <= 1
 sol(3109): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3110): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,2,13) <= 1
 sol(3111): lambda(1,1,14) + lambda(1,3,8) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(3112): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3113): lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(3114): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(3115): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3116): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(6,2,11) + lambda(3,3,9) <= 1
 sol(3117): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3118): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3119): lambda(1,1,8) + lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(6,2,11) <= 1
 sol(3120): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) + lambda(3,3,9) <= 1
 sol(3121): lambda(6,4,11) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3122): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3123): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,15) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3124): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(1,2,15) <= 1
 sol(3125): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3126): lambda(4,1,13) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3127): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3128): lambda(3,3,13) + lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(5,2,13) <= 1
 sol(3129): lambda(2,2,8) + lambda(6,1,15) + lambda(1,3,8) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,3,10) <= 1
 sol(3130): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3131): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3132): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3133): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3134): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3135): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(3136): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(1,2,15) <= 1
 sol(3137): lambda(1,1,14) + lambda(6,6,8) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3138): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(3,6,9)
   + lambda(4,2,8) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(3139): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3140): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3141): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,13) + lambda(4,1,8) + lambda(4,2,13)
   + lambda(1,3,10) <= 1
 sol(3142): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3143): lambda(4,1,13) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(6,4,14) + lambda(5,3,15)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3144): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3145): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3146): lambda(1,1,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(2,3,8) + lambda(2,3,13) + lambda(4,2,13) <= 1
 sol(3147): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3148): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(5,2,13)
   + lambda(1,3,10) <= 1
 sol(3149): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(4,2,8) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3150): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(3151): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3152): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,13) + lambda(6,2,11)
   <= 1
 sol(3153): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(3154): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3155): lambda(6,6,9) + lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3156): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,4,14)
   + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(3157): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3158): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,3,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(6,2,11) <= 1
 sol(3159): lambda(2,2,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,3,10) <= 1
 sol(3160): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3161): lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) <= 1
 sol(3162): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3163): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,3,11) + lambda(3,6,9) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(4,2,13) <= 1
 sol(3164): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3165): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(2,3,11) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(3166): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3167): lambda(3,3,11) + lambda(1,1,14) + lambda(6,6,8) + lambda(1,3,8)
   + lambda(2,5,13) + lambda(4,1,8) + lambda(4,2,13) + lambda(1,3,10) <= 1
 sol(3168): lambda(1,1,8) + lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3169): lambda(6,4,11) + lambda(3,3,13) + lambda(4,1,13)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3170): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,6,15)
   + lambda(3,4,11) + lambda(2,3,8) <= 1
 sol(3171): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(4,2,13) <= 1
 sol(3172): lambda(6,1,15) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(5,2,13) <= 1
 sol(3173): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3174): lambda(4,1,13) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3175): lambda(1,1,14) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(6,6,14) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3176): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3177): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3178): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(2,3,8) <= 1
 sol(3179): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(4,1,8) + lambda(1,2,15)
   + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(3180): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(3,4,11)
   + lambda(4,3,14) + lambda(2,3,8) <= 1
 sol(3181): lambda(2,2,8) + lambda(1,1,8) + lambda(4,1,13) + lambda(3,6,15)
   + lambda(5,3,15) + lambda(3,4,11) <= 1
 sol(3182): lambda(6,1,15) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3183): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3184): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3185): lambda(6,4,11) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,11) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3186): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,11) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3187): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(5,3,13)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3188): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(1,2,15) <= 1
 sol(3189): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(1,3,10) + lambda(6,2,11) <= 1
 sol(3190): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(2,3,8)
   + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(3191): lambda(6,4,11) + lambda(3,3,13) + lambda(1,1,14)
   + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3192): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3193): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(6,4,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3194): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,13)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3195): lambda(6,4,11) + lambda(6,1,15) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,11) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(4,1,8) + lambda(1,2,15) <= 1
 sol(3196): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3197): lambda(1,1,8) + lambda(6,4,11) + lambda(4,1,13) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(1,2,15)
   <= 1
 sol(3198): lambda(4,1,13) + lambda(6,6,9) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3199): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(4,2,13) <= 1
 sol(3200): lambda(3,3,13) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(6,4,14) + lambda(2,3,8) + lambda(4,1,8)
   + lambda(6,2,11) <= 1
 sol(3201): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(5,2,13) <= 1
 sol(3202): lambda(6,4,11) + lambda(1,1,9) + lambda(2,5,13)
   + lambda(3,4,14) + lambda(4,2,8) + lambda(6,6,14) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3203): lambda(4,1,13) + lambda(6,6,8) + lambda(1,3,8) + lambda(2,3,11)
   + lambda(3,4,14) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10) <= 1
 sol(3204): lambda(1,1,8) + lambda(6,4,11) + lambda(3,3,13)
   + lambda(4,1,13) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3205): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3206): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(5,3,13) + lambda(3,6,15) <= 1
 sol(3207): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(4,2,8) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3208): lambda(6,4,11) + lambda(1,1,9) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3209): lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,10) + lambda(6,4,14)
   + lambda(4,3,14) + lambda(5,2,13) <= 1
 sol(3210): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(4,3,14)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3211): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(2,3,8) + lambda(1,2,15) <= 1
 sol(3212): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3213): lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8) + lambda(2,3,11)
   + lambda(3,6,9) + lambda(4,1,8) + lambda(1,2,15) + lambda(1,3,10)
   + lambda(6,2,11) <= 1
 sol(3214): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3215): lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,4,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(1,2,15) + lambda(3,3,9) <= 1
 sol(3216): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3217): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,9) + lambda(4,3,14) <= 1
 sol(3218): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,9) + lambda(5,3,15) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3219): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3220): lambda(6,4,11) + lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,8)
   + lambda(2,3,11) + lambda(3,6,10) + lambda(6,4,14) + lambda(4,1,8)
   + lambda(5,2,13) + lambda(1,3,10) <= 1
 sol(3221): lambda(2,2,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(1,3,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(5,3,15) + lambda(4,1,8)
   + lambda(3,3,9) <= 1
 sol(3222): lambda(6,4,11) + lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(4,2,13) <= 1
 sol(3223): lambda(1,1,8) + lambda(6,4,11) + lambda(2,5,13) + lambda(4,2,8)
   + lambda(3,6,15) + lambda(6,4,14) + lambda(4,3,14) + lambda(1,2,15)
   + lambda(3,3,9) <= 1
 sol(3224): lambda(4,1,13) + lambda(1,3,8) + lambda(6,2,13)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(3,6,15) + lambda(5,3,15)
   + lambda(3,4,11) + lambda(1,3,10) <= 1
 sol(3225): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3226): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,9) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3227): lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3228): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(5,3,13) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3229): lambda(2,2,8) + lambda(6,4,11) + lambda(6,1,15) + lambda(1,3,8)
   + lambda(1,1,9) + lambda(3,6,15) + lambda(5,3,15) + lambda(3,3,9) <= 1
 sol(3230): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(6,4,14) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(6,2,11) <= 1
 sol(3231): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,10) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3232): lambda(6,4,11) + lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(2,5,13) + lambda(3,6,10) + lambda(5,3,15) + lambda(2,3,8)
   + lambda(1,2,15) <= 1
 sol(3233): lambda(1,1,8) + lambda(1,1,14) + lambda(6,2,8) + lambda(3,6,10)
   + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13) + lambda(6,2,11) <= 1
 sol(3234): lambda(1,1,8) + lambda(6,4,11) + lambda(1,1,14) + lambda(6,2,8)
   + lambda(3,4,14) + lambda(6,6,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(5,2,13) <= 1
 sol(3235): lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13) + lambda(6,6,14)
   + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8) + lambda(1,2,15)
   + lambda(6,2,11) <= 1
 sol(3236): lambda(1,1,8) + lambda(6,6,9) + lambda(1,1,14) + lambda(4,2,8)
   + lambda(3,4,11) + lambda(2,3,8) + lambda(2,3,13) + lambda(5,2,13) <= 1
 sol(3237): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3238): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(5,3,13) + lambda(3,6,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3239): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(5,3,15) + lambda(3,4,11) + lambda(2,3,8)
   + lambda(4,2,13) <= 1
 sol(3240): lambda(6,6,9) + lambda(1,1,14) + lambda(1,3,8) + lambda(1,1,9)
   + lambda(2,5,13) + lambda(4,2,8) + lambda(5,3,13) + lambda(3,4,11)
   + lambda(4,2,13) <= 1
 sol(3241): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(6,6,14) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3242): lambda(1,1,14) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
 sol(3243): lambda(6,4,11) + lambda(4,1,13) + lambda(1,1,9) + lambda(6,2,8)
   + lambda(3,6,9) + lambda(6,4,14) + lambda(2,3,8) + lambda(2,3,13)
   + lambda(1,2,15) <= 1
 sol(3244): lambda(6,1,15) + lambda(1,1,9) + lambda(6,2,8) + lambda(2,5,13)
   + lambda(3,6,15) + lambda(3,4,11) + lambda(4,3,14) + lambda(2,3,8)
   + lambda(5,2,13) <= 1
Bounds
 lambda(2,2,8) <= 1
 lambda(3,3,11) <= 1
 lambda(1,1,8) <= 1
 lambda(6,4,11) <= 1
 lambda(6,1,15) <= 1
 lambda(3,3,13) <= 1
 lambda(4,1,13) <= 1
 lambda(6,6,9) <= 1
 lambda(1,1,14) <= 1
 lambda(6,6,8) <= 1
 lambda(1,3,8) <= 1
 lambda(1,1,9) <= 1
 lambda(6,2,8) <= 1
 lambda(6,2,13) <= 1
 lambda(2,5,13) <= 1
 lambda(2,3,11) <= 1
 lambda(3,4,14) <= 1
 lambda(3,6,9) <= 1
 lambda(4,2,8) <= 1
 lambda(3,6,11) <= 1
 lambda(3,6,10) <= 1
 lambda(5,3,13) <= 1
 lambda(3,6,15) <= 1
 lambda(6,6,14) <= 1
 lambda(6,4,14) <= 1
 lambda(5,3,15) <= 1
 lambda(3,4,11) <= 1
 lambda(4,3,14) <= 1
 lambda(2,3,8) <= 1
 lambda(2,3,13) <= 1
 lambda(4,1,8) <= 1
 lambda(5,2,13) <= 1
 lambda(4,2,13) <= 1
 lambda(1,2,15) <= 1
 lambda(1,3,10) <= 1
 lambda(6,2,11) <= 1
 lambda(3,3,9) <= 1
End
