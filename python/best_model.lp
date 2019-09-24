\Problem name: 

Minimize
OBJROW: - xij(1,0) - xij(2,0) - xij(0,1) - xij(2,1) - xij(0,2) - xij(1,2)
Subject To
flow(0):  xij(1,0) + xij(2,0) - xij(0,1) - xij(0,2) = 0
flow(1):  - xij(1,0) + xij(0,1) + xij(2,1) - xij(1,2) = 0
flow(2):  - xij(2,0) - xij(2,1) + xij(0,2) + xij(1,2) = 0
triangle_cut:  - xij(1,0) - xij(2,0) -0.821041921 xij(2,1) -0.478974359 xij(0,2) -0.988481888 xij(1,2) >= -1
Bounds
 0 <= xij(0,0) <= 0
 1 <= xij(1,0) <= 1
 0 <= xij(2,0) <= 0
 1 <= xij(0,1) <= 1
 0 <= xij(1,1) <= 0
 0 <= xij(2,1) <= 0
 0 <= xij(0,2) <= 0
 0 <= xij(1,2) <= 0
 0 <= xij(2,2) <= 0
Integers
xij(0,0) xij(1,0) xij(2,0) xij(0,1) xij(1,1) xij(2,1) xij(0,2) xij(1,2) xij(2,2) 
End
