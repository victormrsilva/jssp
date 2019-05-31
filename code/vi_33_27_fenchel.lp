\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  - 0.2 lambda(2,1,26) - 0.4 lambda(2,1,24) - 0.2 lambda(3,1,21)
   - 0.2 lambda(2,3,25) - 0.4 lambda(2,3,21) - 0.2 lambda(2,2,21)
   - 0.4 lambda(1,1,21) - 0.4 lambda(1,3,19) - 0.4 lambda(3,1,18)
Subject To
 sol(0): lambda(2,1,26) <= 1
 sol(1): lambda(2,1,24) <= 1
 sol(2): lambda(3,1,21) + lambda(2,3,25) <= 1
 sol(3): lambda(2,1,24) + lambda(2,3,21) <= 1
 sol(4): lambda(2,1,26) + lambda(2,3,21) <= 1
 sol(5): lambda(2,1,26) + lambda(2,2,21) <= 1
 sol(6): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) <= 1
 sol(7): lambda(2,1,26) + lambda(1,1,21) <= 1
 sol(8): lambda(2,1,26) + lambda(2,3,25) + lambda(1,1,21) <= 1
 sol(9): lambda(2,1,26) + lambda(2,2,21) + lambda(1,1,21) <= 1
 sol(10): lambda(2,1,26) + lambda(2,3,25) + lambda(1,3,19) <= 1
 sol(11): lambda(2,1,24) + lambda(2,3,21) + lambda(1,3,19) <= 1
 sol(12): lambda(2,1,26) + lambda(2,3,25) + lambda(1,3,19) + lambda(3,1,18)
   <= 1
 sol(13): lambda(2,1,26) + lambda(2,3,21) + lambda(1,1,21) <= 1
 sol(14): lambda(2,1,26) + lambda(1,3,19) <= 1
 sol(15): lambda(2,1,24) + lambda(3,1,18) <= 1
 sol(16): lambda(2,1,26) + lambda(1,3,19) + lambda(3,1,18) <= 1
 sol(17): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) + lambda(1,3,19)
   <= 1
 sol(18): lambda(2,1,24) + lambda(1,3,19) + lambda(3,1,18) <= 1
 sol(19): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) + lambda(1,1,21)
   <= 1
 sol(20): lambda(2,1,24) + lambda(1,3,19) <= 1
 sol(21): lambda(2,1,26) + lambda(2,2,21) + lambda(3,1,18) <= 1
 sol(22): lambda(2,1,24) + lambda(2,3,21) + lambda(3,1,18) <= 1
 sol(23): lambda(2,1,26) + lambda(2,2,21) + lambda(1,1,21) + lambda(1,3,19)
   <= 1
 sol(24): lambda(2,1,26) + lambda(2,3,21) + lambda(1,1,21) + lambda(1,3,19)
   <= 1
 sol(25): lambda(2,1,24) + lambda(2,3,21) + lambda(1,3,19) + lambda(3,1,18)
   <= 1
 sol(26): lambda(2,1,26) + lambda(2,3,21) + lambda(1,3,19) <= 1
 sol(27): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) + lambda(3,1,18)
   <= 1
 sol(28): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) + lambda(1,3,19)
   + lambda(3,1,18) <= 1
 sol(29): lambda(3,1,21) + lambda(2,3,21) <= 1
 sol(30): lambda(3,1,21) + lambda(2,3,25) + lambda(1,3,19) <= 1
 sol(31): lambda(2,1,26) + lambda(1,1,21) + lambda(1,3,19) <= 1
 sol(32): lambda(3,1,21) + lambda(2,3,25) + lambda(2,2,21) <= 1
 sol(33): lambda(2,1,26) + lambda(2,2,21) + lambda(1,3,19) <= 1
 sol(34): lambda(2,1,26) + lambda(2,2,21) + lambda(1,3,19) + lambda(3,1,18)
   <= 1
 sol(35): lambda(2,1,26) + lambda(2,3,25) + lambda(2,2,21) + lambda(1,1,21)
   + lambda(1,3,19) <= 1
 sol(36): lambda(2,1,26) + lambda(2,3,21) + lambda(3,1,18) <= 1
 sol(37): lambda(3,1,21) + lambda(2,3,21) + lambda(1,3,19) <= 1
 sol(38): lambda(2,1,26) + lambda(2,3,25) + lambda(3,1,18) <= 1
 sol(39): lambda(2,1,26) + lambda(2,3,21) + lambda(1,3,19) + lambda(3,1,18)
   <= 1
 sol(40): lambda(2,1,26) + lambda(3,1,18) <= 1
 sol(41): lambda(2,1,26) + lambda(2,3,25) <= 1
 sol(42): lambda(2,1,26) + lambda(2,3,25) + lambda(1,1,21) + lambda(1,3,19)
   <= 1
 sol(43): lambda(3,1,21) + lambda(2,3,25) + lambda(2,2,21) + lambda(1,3,19)
   <= 1
Bounds
 lambda(2,1,26) <= 1
 lambda(2,1,24) <= 1
 lambda(3,1,21) <= 1
 lambda(2,3,25) <= 1
 lambda(2,3,21) <= 1
 lambda(2,2,21) <= 1
 lambda(1,1,21) <= 1
 lambda(1,3,19) <= 1
 lambda(3,1,18) <= 1
End
