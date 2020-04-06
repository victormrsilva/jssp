\Problem name: 

Minimize
OBJROW: - xij(0,1) - xij(0,2) - xij(1,0) - xij(1,2) - xij(2,0) - xij(2,1)
Subject To
minimum_3:  xij(0,1) + xij(0,2) + xij(1,0) + xij(1,2) + xij(2,0) + xij(2,1) >= 3
flow(0):  xij(0,1) + xij(0,2) - xij(1,0) - xij(2,0) = 0
flow(1):  - xij(0,1) + xij(1,0) + xij(1,2) - xij(2,1) = 0
flow(2):  - xij(0,2) - xij(1,2) + xij(2,0) + xij(2,1) = 0
triangle_cut:  - xij(0,2) -0.25000 xij(1,0) -0.75000 xij(1,2) - xij(2,0) - xij(2,1) >= -1
Bounds
 0 <= xij(0,1) <= 1
 0 <= xij(0,2) <= 1
 0 <= xij(1,0) <= 1
 0 <= xij(1,2) <= 1
 0 <= xij(2,0) <= 1
 0 <= xij(2,1) <= 1
Integers
xij(0,1) xij(0,2) xij(1,0) xij(1,2) xij(2,0) xij(2,1) 
End
