\ Model compact
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C
Subject To
 ord(0,1): x(j0,m0) - x(j0,m2) >= 1
 ord(0,2): - x(j0,m0) + x(j0,m1) >= 3
 ord(0,3): - x(j0,m1) + x(j0,m3) >= 6
 ord(0,4): - x(j0,m3) + x(j0,m5) >= 7
 ord(0,5): x(j0,m4) - x(j0,m5) >= 3
 ord(1,1): - x(j1,m1) + x(j1,m2) >= 8
 ord(1,2): - x(j1,m2) + x(j1,m4) >= 5
 ord(1,3): - x(j1,m4) + x(j1,m5) >= 10
 ord(1,4): x(j1,m0) - x(j1,m5) >= 10
 ord(1,5): - x(j1,m0) + x(j1,m3) >= 10
 ord(2,1): - x(j2,m2) + x(j2,m3) >= 5
 ord(2,2): - x(j2,m3) + x(j2,m5) >= 4
 ord(2,3): x(j2,m0) - x(j2,m5) >= 8
 ord(2,4): - x(j2,m0) + x(j2,m1) >= 9
 ord(2,5): - x(j2,m1) + x(j2,m4) >= 1
 ord(3,1): x(j3,m0) - x(j3,m1) >= 5
 ord(3,2): - x(j3,m0) + x(j3,m2) >= 5
 ord(3,3): - x(j3,m2) + x(j3,m3) >= 5
 ord(3,4): - x(j3,m3) + x(j3,m4) >= 3
 ord(3,5): - x(j3,m4) + x(j3,m5) >= 8
 ord(4,1): x(j4,m1) - x(j4,m2) >= 9
 ord(4,2): - x(j4,m1) + x(j4,m4) >= 3
 ord(4,3): - x(j4,m4) + x(j4,m5) >= 5
 ord(4,4): x(j4,m0) - x(j4,m5) >= 4
 ord(4,5): - x(j4,m0) + x(j4,m3) >= 3
 ord(5,1): - x(j5,m1) + x(j5,m3) >= 3
 ord(5,2): - x(j5,m3) + x(j5,m5) >= 3
 ord(5,3): x(j5,m0) - x(j5,m5) >= 9
 ord(5,4): - x(j5,m0) + x(j5,m4) >= 10
 ord(5,5): x(j5,m2) - x(j5,m4) >= 4
 psy(0,1,0): x(j0,m0) - 197 y(j0,k1,m0) - x(j1,m0) >= -187
 psy(0,2,0): x(j0,m0) - 197 y(j0,k2,m0) - x(j2,m0) >= -188
 psy(0,3,0): x(j0,m0) - 197 y(j0,k3,m0) - x(j3,m0) >= -192
 psy(0,4,0): x(j0,m0) - 197 y(j0,k4,m0) - x(j4,m0) >= -194
 psy(0,5,0): x(j0,m0) - 197 y(j0,k5,m0) - x(j5,m0) >= -187
 sum(0,0): y(j0,k1,m0) + y(j0,k2,m0) + y(j0,k3,m0) + y(j0,k4,m0)
   + y(j0,k5,m0) + y(j0,k6,m0) = 1
 psy(1,0,0): - x(j0,m0) - 197 y(j1,k0,m0) + x(j1,m0) >= -194
 psy(1,2,0): x(j1,m0) - 197 y(j1,k2,m0) - x(j2,m0) >= -188
 psy(1,3,0): x(j1,m0) - 197 y(j1,k3,m0) - x(j3,m0) >= -192
 psy(1,4,0): x(j1,m0) - 197 y(j1,k4,m0) - x(j4,m0) >= -194
 psy(1,5,0): x(j1,m0) - 197 y(j1,k5,m0) - x(j5,m0) >= -187
 sum(1,0): y(j1,k0,m0) + y(j1,k2,m0) + y(j1,k3,m0) + y(j1,k4,m0)
   + y(j1,k5,m0) + y(j1,k6,m0) = 1
 psy(2,0,0): - x(j0,m0) - 197 y(j2,k0,m0) + x(j2,m0) >= -194
 psy(2,1,0): - x(j1,m0) - 197 y(j2,k1,m0) + x(j2,m0) >= -187
 psy(2,3,0): x(j2,m0) - 197 y(j2,k3,m0) - x(j3,m0) >= -192
 psy(2,4,0): x(j2,m0) - 197 y(j2,k4,m0) - x(j4,m0) >= -194
 psy(2,5,0): x(j2,m0) - 197 y(j2,k5,m0) - x(j5,m0) >= -187
 sum(2,0): y(j2,k0,m0) + y(j2,k1,m0) + y(j2,k3,m0) + y(j2,k4,m0)
   + y(j2,k5,m0) + y(j2,k6,m0) = 1
 psy(3,0,0): - x(j0,m0) - 197 y(j3,k0,m0) + x(j3,m0) >= -194
 psy(3,1,0): - x(j1,m0) - 197 y(j3,k1,m0) + x(j3,m0) >= -187
 psy(3,2,0): - x(j2,m0) - 197 y(j3,k2,m0) + x(j3,m0) >= -188
 psy(3,4,0): x(j3,m0) - 197 y(j3,k4,m0) - x(j4,m0) >= -194
 psy(3,5,0): x(j3,m0) - 197 y(j3,k5,m0) - x(j5,m0) >= -187
 sum(3,0): y(j3,k0,m0) + y(j3,k1,m0) + y(j3,k2,m0) + y(j3,k4,m0)
   + y(j3,k5,m0) + y(j3,k6,m0) = 1
 psy(4,0,0): - x(j0,m0) - 197 y(j4,k0,m0) + x(j4,m0) >= -194
 psy(4,1,0): - x(j1,m0) - 197 y(j4,k1,m0) + x(j4,m0) >= -187
 psy(4,2,0): - x(j2,m0) - 197 y(j4,k2,m0) + x(j4,m0) >= -188
 psy(4,3,0): - x(j3,m0) - 197 y(j4,k3,m0) + x(j4,m0) >= -192
 psy(4,5,0): x(j4,m0) - 197 y(j4,k5,m0) - x(j5,m0) >= -187
 sum(4,0): y(j4,k0,m0) + y(j4,k1,m0) + y(j4,k2,m0) + y(j4,k3,m0)
   + y(j4,k5,m0) + y(j4,k6,m0) = 1
 psy(5,0,0): - x(j0,m0) - 197 y(j5,k0,m0) + x(j5,m0) >= -194
 psy(5,1,0): - x(j1,m0) - 197 y(j5,k1,m0) + x(j5,m0) >= -187
 psy(5,2,0): - x(j2,m0) - 197 y(j5,k2,m0) + x(j5,m0) >= -188
 psy(5,3,0): - x(j3,m0) - 197 y(j5,k3,m0) + x(j5,m0) >= -192
 psy(5,4,0): - x(j4,m0) - 197 y(j5,k4,m0) + x(j5,m0) >= -194
 sum(5,0): y(j5,k0,m0) + y(j5,k1,m0) + y(j5,k2,m0) + y(j5,k3,m0)
   + y(j5,k4,m0) + y(j5,k6,m0) = 1
 sum2(0,0): y(j1,k0,m0) + y(j2,k0,m0) + y(j3,k0,m0) + y(j4,k0,m0)
   + y(j5,k0,m0) <= 1
 sum2(1,0): y(j0,k1,m0) + y(j2,k1,m0) + y(j3,k1,m0) + y(j4,k1,m0)
   + y(j5,k1,m0) <= 1
 sum2(2,0): y(j0,k2,m0) + y(j1,k2,m0) + y(j3,k2,m0) + y(j4,k2,m0)
   + y(j5,k2,m0) <= 1
 sum2(3,0): y(j0,k3,m0) + y(j1,k3,m0) + y(j2,k3,m0) + y(j4,k3,m0)
   + y(j5,k3,m0) <= 1
 sum2(4,0): y(j0,k4,m0) + y(j1,k4,m0) + y(j2,k4,m0) + y(j3,k4,m0)
   + y(j5,k4,m0) <= 1
 sum2(5,0): y(j0,k5,m0) + y(j1,k5,m0) + y(j2,k5,m0) + y(j3,k5,m0)
   + y(j4,k5,m0) <= 1
 sum2(6,0): y(j0,k6,m0) + y(j1,k6,m0) + y(j2,k6,m0) + y(j3,k6,m0)
   + y(j4,k6,m0) + y(j5,k6,m0) <= 1
 psy(0,1,1): x(j0,m1) - 197 y(j0,k1,m1) - x(j1,m1) >= -189
 psy(0,2,1): x(j0,m1) - 197 y(j0,k2,m1) - x(j2,m1) >= -196
 psy(0,3,1): x(j0,m1) - 197 y(j0,k3,m1) - x(j3,m1) >= -192
 psy(0,4,1): x(j0,m1) - 197 y(j0,k4,m1) - x(j4,m1) >= -194
 psy(0,5,1): x(j0,m1) - 197 y(j0,k5,m1) - x(j5,m1) >= -194
 sum(0,1): y(j0,k1,m1) + y(j0,k2,m1) + y(j0,k3,m1) + y(j0,k4,m1)
   + y(j0,k5,m1) + y(j0,k6,m1) = 1
 psy(1,0,1): - x(j0,m1) - 197 y(j1,k0,m1) + x(j1,m1) >= -191
 psy(1,2,1): x(j1,m1) - 197 y(j1,k2,m1) - x(j2,m1) >= -196
 psy(1,3,1): x(j1,m1) - 197 y(j1,k3,m1) - x(j3,m1) >= -192
 psy(1,4,1): x(j1,m1) - 197 y(j1,k4,m1) - x(j4,m1) >= -194
 psy(1,5,1): x(j1,m1) - 197 y(j1,k5,m1) - x(j5,m1) >= -194
 sum(1,1): y(j1,k0,m1) + y(j1,k2,m1) + y(j1,k3,m1) + y(j1,k4,m1)
   + y(j1,k5,m1) + y(j1,k6,m1) = 1
 psy(2,0,1): - x(j0,m1) - 197 y(j2,k0,m1) + x(j2,m1) >= -191
 psy(2,1,1): - x(j1,m1) - 197 y(j2,k1,m1) + x(j2,m1) >= -189
 psy(2,3,1): x(j2,m1) - 197 y(j2,k3,m1) - x(j3,m1) >= -192
 psy(2,4,1): x(j2,m1) - 197 y(j2,k4,m1) - x(j4,m1) >= -194
 psy(2,5,1): x(j2,m1) - 197 y(j2,k5,m1) - x(j5,m1) >= -194
 sum(2,1): y(j2,k0,m1) + y(j2,k1,m1) + y(j2,k3,m1) + y(j2,k4,m1)
   + y(j2,k5,m1) + y(j2,k6,m1) = 1
 psy(3,0,1): - x(j0,m1) - 197 y(j3,k0,m1) + x(j3,m1) >= -191
 psy(3,1,1): - x(j1,m1) - 197 y(j3,k1,m1) + x(j3,m1) >= -189
 psy(3,2,1): - x(j2,m1) - 197 y(j3,k2,m1) + x(j3,m1) >= -196
 psy(3,4,1): x(j3,m1) - 197 y(j3,k4,m1) - x(j4,m1) >= -194
 psy(3,5,1): x(j3,m1) - 197 y(j3,k5,m1) - x(j5,m1) >= -194
 sum(3,1): y(j3,k0,m1) + y(j3,k1,m1) + y(j3,k2,m1) + y(j3,k4,m1)
   + y(j3,k5,m1) + y(j3,k6,m1) = 1
 psy(4,0,1): - x(j0,m1) - 197 y(j4,k0,m1) + x(j4,m1) >= -191
 psy(4,1,1): - x(j1,m1) - 197 y(j4,k1,m1) + x(j4,m1) >= -189
 psy(4,2,1): - x(j2,m1) - 197 y(j4,k2,m1) + x(j4,m1) >= -196
 psy(4,3,1): - x(j3,m1) - 197 y(j4,k3,m1) + x(j4,m1) >= -192
 psy(4,5,1): x(j4,m1) - 197 y(j4,k5,m1) - x(j5,m1) >= -194
 sum(4,1): y(j4,k0,m1) + y(j4,k1,m1) + y(j4,k2,m1) + y(j4,k3,m1)
   + y(j4,k5,m1) + y(j4,k6,m1) = 1
 psy(5,0,1): - x(j0,m1) - 197 y(j5,k0,m1) + x(j5,m1) >= -191
 psy(5,1,1): - x(j1,m1) - 197 y(j5,k1,m1) + x(j5,m1) >= -189
 psy(5,2,1): - x(j2,m1) - 197 y(j5,k2,m1) + x(j5,m1) >= -196
 psy(5,3,1): - x(j3,m1) - 197 y(j5,k3,m1) + x(j5,m1) >= -192
 psy(5,4,1): - x(j4,m1) - 197 y(j5,k4,m1) + x(j5,m1) >= -194
 sum(5,1): y(j5,k0,m1) + y(j5,k1,m1) + y(j5,k2,m1) + y(j5,k3,m1)
   + y(j5,k4,m1) + y(j5,k6,m1) = 1
 sum2(0,1): y(j1,k0,m1) + y(j2,k0,m1) + y(j3,k0,m1) + y(j4,k0,m1)
   + y(j5,k0,m1) <= 1
 sum2(1,1): y(j0,k1,m1) + y(j2,k1,m1) + y(j3,k1,m1) + y(j4,k1,m1)
   + y(j5,k1,m1) <= 1
 sum2(2,1): y(j0,k2,m1) + y(j1,k2,m1) + y(j3,k2,m1) + y(j4,k2,m1)
   + y(j5,k2,m1) <= 1
 sum2(3,1): y(j0,k3,m1) + y(j1,k3,m1) + y(j2,k3,m1) + y(j4,k3,m1)
   + y(j5,k3,m1) <= 1
 sum2(4,1): y(j0,k4,m1) + y(j1,k4,m1) + y(j2,k4,m1) + y(j3,k4,m1)
   + y(j5,k4,m1) <= 1
 sum2(5,1): y(j0,k5,m1) + y(j1,k5,m1) + y(j2,k5,m1) + y(j3,k5,m1)
   + y(j4,k5,m1) <= 1
 sum2(6,1): y(j0,k6,m1) + y(j1,k6,m1) + y(j2,k6,m1) + y(j3,k6,m1)
   + y(j4,k6,m1) + y(j5,k6,m1) <= 1
 psy(0,1,2): x(j0,m2) - 197 y(j0,k1,m2) - x(j1,m2) >= -192
 psy(0,2,2): x(j0,m2) - 197 y(j0,k2,m2) - x(j2,m2) >= -192
 psy(0,3,2): x(j0,m2) - 197 y(j0,k3,m2) - x(j3,m2) >= -192
 psy(0,4,2): x(j0,m2) - 197 y(j0,k4,m2) - x(j4,m2) >= -188
 psy(0,5,2): x(j0,m2) - 197 y(j0,k5,m2) - x(j5,m2) >= -196
 sum(0,2): y(j0,k1,m2) + y(j0,k2,m2) + y(j0,k3,m2) + y(j0,k4,m2)
   + y(j0,k5,m2) + y(j0,k6,m2) = 1
 psy(1,0,2): - x(j0,m2) - 197 y(j1,k0,m2) + x(j1,m2) >= -196
 psy(1,2,2): x(j1,m2) - 197 y(j1,k2,m2) - x(j2,m2) >= -192
 psy(1,3,2): x(j1,m2) - 197 y(j1,k3,m2) - x(j3,m2) >= -192
 psy(1,4,2): x(j1,m2) - 197 y(j1,k4,m2) - x(j4,m2) >= -188
 psy(1,5,2): x(j1,m2) - 197 y(j1,k5,m2) - x(j5,m2) >= -196
 sum(1,2): y(j1,k0,m2) + y(j1,k2,m2) + y(j1,k3,m2) + y(j1,k4,m2)
   + y(j1,k5,m2) + y(j1,k6,m2) = 1
 psy(2,0,2): - x(j0,m2) - 197 y(j2,k0,m2) + x(j2,m2) >= -196
 psy(2,1,2): - x(j1,m2) - 197 y(j2,k1,m2) + x(j2,m2) >= -192
 psy(2,3,2): x(j2,m2) - 197 y(j2,k3,m2) - x(j3,m2) >= -192
 psy(2,4,2): x(j2,m2) - 197 y(j2,k4,m2) - x(j4,m2) >= -188
 psy(2,5,2): x(j2,m2) - 197 y(j2,k5,m2) - x(j5,m2) >= -196
 sum(2,2): y(j2,k0,m2) + y(j2,k1,m2) + y(j2,k3,m2) + y(j2,k4,m2)
   + y(j2,k5,m2) + y(j2,k6,m2) = 1
 psy(3,0,2): - x(j0,m2) - 197 y(j3,k0,m2) + x(j3,m2) >= -196
 psy(3,1,2): - x(j1,m2) - 197 y(j3,k1,m2) + x(j3,m2) >= -192
 psy(3,2,2): - x(j2,m2) - 197 y(j3,k2,m2) + x(j3,m2) >= -192
 psy(3,4,2): x(j3,m2) - 197 y(j3,k4,m2) - x(j4,m2) >= -188
 psy(3,5,2): x(j3,m2) - 197 y(j3,k5,m2) - x(j5,m2) >= -196
 sum(3,2): y(j3,k0,m2) + y(j3,k1,m2) + y(j3,k2,m2) + y(j3,k4,m2)
   + y(j3,k5,m2) + y(j3,k6,m2) = 1
 psy(4,0,2): - x(j0,m2) - 197 y(j4,k0,m2) + x(j4,m2) >= -196
 psy(4,1,2): - x(j1,m2) - 197 y(j4,k1,m2) + x(j4,m2) >= -192
 psy(4,2,2): - x(j2,m2) - 197 y(j4,k2,m2) + x(j4,m2) >= -192
 psy(4,3,2): - x(j3,m2) - 197 y(j4,k3,m2) + x(j4,m2) >= -192
 psy(4,5,2): x(j4,m2) - 197 y(j4,k5,m2) - x(j5,m2) >= -196
 sum(4,2): y(j4,k0,m2) + y(j4,k1,m2) + y(j4,k2,m2) + y(j4,k3,m2)
   + y(j4,k5,m2) + y(j4,k6,m2) = 1
 psy(5,0,2): - x(j0,m2) - 197 y(j5,k0,m2) + x(j5,m2) >= -196
 psy(5,1,2): - x(j1,m2) - 197 y(j5,k1,m2) + x(j5,m2) >= -192
 psy(5,2,2): - x(j2,m2) - 197 y(j5,k2,m2) + x(j5,m2) >= -192
 psy(5,3,2): - x(j3,m2) - 197 y(j5,k3,m2) + x(j5,m2) >= -192
 psy(5,4,2): - x(j4,m2) - 197 y(j5,k4,m2) + x(j5,m2) >= -188
 sum(5,2): y(j5,k0,m2) + y(j5,k1,m2) + y(j5,k2,m2) + y(j5,k3,m2)
   + y(j5,k4,m2) + y(j5,k6,m2) = 1
 sum2(0,2): y(j1,k0,m2) + y(j2,k0,m2) + y(j3,k0,m2) + y(j4,k0,m2)
   + y(j5,k0,m2) <= 1
 sum2(1,2): y(j0,k1,m2) + y(j2,k1,m2) + y(j3,k1,m2) + y(j4,k1,m2)
   + y(j5,k1,m2) <= 1
 sum2(2,2): y(j0,k2,m2) + y(j1,k2,m2) + y(j3,k2,m2) + y(j4,k2,m2)
   + y(j5,k2,m2) <= 1
 sum2(3,2): y(j0,k3,m2) + y(j1,k3,m2) + y(j2,k3,m2) + y(j4,k3,m2)
   + y(j5,k3,m2) <= 1
 sum2(4,2): y(j0,k4,m2) + y(j1,k4,m2) + y(j2,k4,m2) + y(j3,k4,m2)
   + y(j5,k4,m2) <= 1
 sum2(5,2): y(j0,k5,m2) + y(j1,k5,m2) + y(j2,k5,m2) + y(j3,k5,m2)
   + y(j4,k5,m2) <= 1
 sum2(6,2): y(j0,k6,m2) + y(j1,k6,m2) + y(j2,k6,m2) + y(j3,k6,m2)
   + y(j4,k6,m2) + y(j5,k6,m2) <= 1
 psy(0,1,3): x(j0,m3) - 197 y(j0,k1,m3) - x(j1,m3) >= -193
 psy(0,2,3): x(j0,m3) - 197 y(j0,k2,m3) - x(j2,m3) >= -193
 psy(0,3,3): x(j0,m3) - 197 y(j0,k3,m3) - x(j3,m3) >= -194
 psy(0,4,3): x(j0,m3) - 197 y(j0,k4,m3) - x(j4,m3) >= -196
 psy(0,5,3): x(j0,m3) - 197 y(j0,k5,m3) - x(j5,m3) >= -194
 sum(0,3): y(j0,k1,m3) + y(j0,k2,m3) + y(j0,k3,m3) + y(j0,k4,m3)
   + y(j0,k5,m3) + y(j0,k6,m3) = 1
 psy(1,0,3): - x(j0,m3) - 197 y(j1,k0,m3) + x(j1,m3) >= -190
 psy(1,2,3): x(j1,m3) - 197 y(j1,k2,m3) - x(j2,m3) >= -193
 psy(1,3,3): x(j1,m3) - 197 y(j1,k3,m3) - x(j3,m3) >= -194
 psy(1,4,3): x(j1,m3) - 197 y(j1,k4,m3) - x(j4,m3) >= -196
 psy(1,5,3): x(j1,m3) - 197 y(j1,k5,m3) - x(j5,m3) >= -194
 sum(1,3): y(j1,k0,m3) + y(j1,k2,m3) + y(j1,k3,m3) + y(j1,k4,m3)
   + y(j1,k5,m3) + y(j1,k6,m3) = 1
 psy(2,0,3): - x(j0,m3) - 197 y(j2,k0,m3) + x(j2,m3) >= -190
 psy(2,1,3): - x(j1,m3) - 197 y(j2,k1,m3) + x(j2,m3) >= -193
 psy(2,3,3): x(j2,m3) - 197 y(j2,k3,m3) - x(j3,m3) >= -194
 psy(2,4,3): x(j2,m3) - 197 y(j2,k4,m3) - x(j4,m3) >= -196
 psy(2,5,3): x(j2,m3) - 197 y(j2,k5,m3) - x(j5,m3) >= -194
 sum(2,3): y(j2,k0,m3) + y(j2,k1,m3) + y(j2,k3,m3) + y(j2,k4,m3)
   + y(j2,k5,m3) + y(j2,k6,m3) = 1
 psy(3,0,3): - x(j0,m3) - 197 y(j3,k0,m3) + x(j3,m3) >= -190
 psy(3,1,3): - x(j1,m3) - 197 y(j3,k1,m3) + x(j3,m3) >= -193
 psy(3,2,3): - x(j2,m3) - 197 y(j3,k2,m3) + x(j3,m3) >= -193
 psy(3,4,3): x(j3,m3) - 197 y(j3,k4,m3) - x(j4,m3) >= -196
 psy(3,5,3): x(j3,m3) - 197 y(j3,k5,m3) - x(j5,m3) >= -194
 sum(3,3): y(j3,k0,m3) + y(j3,k1,m3) + y(j3,k2,m3) + y(j3,k4,m3)
   + y(j3,k5,m3) + y(j3,k6,m3) = 1
 psy(4,0,3): - x(j0,m3) - 197 y(j4,k0,m3) + x(j4,m3) >= -190
 psy(4,1,3): - x(j1,m3) - 197 y(j4,k1,m3) + x(j4,m3) >= -193
 psy(4,2,3): - x(j2,m3) - 197 y(j4,k2,m3) + x(j4,m3) >= -193
 psy(4,3,3): - x(j3,m3) - 197 y(j4,k3,m3) + x(j4,m3) >= -194
 psy(4,5,3): x(j4,m3) - 197 y(j4,k5,m3) - x(j5,m3) >= -194
 sum(4,3): y(j4,k0,m3) + y(j4,k1,m3) + y(j4,k2,m3) + y(j4,k3,m3)
   + y(j4,k5,m3) + y(j4,k6,m3) = 1
 psy(5,0,3): - x(j0,m3) - 197 y(j5,k0,m3) + x(j5,m3) >= -190
 psy(5,1,3): - x(j1,m3) - 197 y(j5,k1,m3) + x(j5,m3) >= -193
 psy(5,2,3): - x(j2,m3) - 197 y(j5,k2,m3) + x(j5,m3) >= -193
 psy(5,3,3): - x(j3,m3) - 197 y(j5,k3,m3) + x(j5,m3) >= -194
 psy(5,4,3): - x(j4,m3) - 197 y(j5,k4,m3) + x(j5,m3) >= -196
 sum(5,3): y(j5,k0,m3) + y(j5,k1,m3) + y(j5,k2,m3) + y(j5,k3,m3)
   + y(j5,k4,m3) + y(j5,k6,m3) = 1
 sum2(0,3): y(j1,k0,m3) + y(j2,k0,m3) + y(j3,k0,m3) + y(j4,k0,m3)
   + y(j5,k0,m3) <= 1
 sum2(1,3): y(j0,k1,m3) + y(j2,k1,m3) + y(j3,k1,m3) + y(j4,k1,m3)
   + y(j5,k1,m3) <= 1
 sum2(2,3): y(j0,k2,m3) + y(j1,k2,m3) + y(j3,k2,m3) + y(j4,k2,m3)
   + y(j5,k2,m3) <= 1
 sum2(3,3): y(j0,k3,m3) + y(j1,k3,m3) + y(j2,k3,m3) + y(j4,k3,m3)
   + y(j5,k3,m3) <= 1
 sum2(4,3): y(j0,k4,m3) + y(j1,k4,m3) + y(j2,k4,m3) + y(j3,k4,m3)
   + y(j5,k4,m3) <= 1
 sum2(5,3): y(j0,k5,m3) + y(j1,k5,m3) + y(j2,k5,m3) + y(j3,k5,m3)
   + y(j4,k5,m3) <= 1
 sum2(6,3): y(j0,k6,m3) + y(j1,k6,m3) + y(j2,k6,m3) + y(j3,k6,m3)
   + y(j4,k6,m3) + y(j5,k6,m3) <= 1
 psy(0,1,4): x(j0,m4) - 197 y(j0,k1,m4) - x(j1,m4) >= -187
 psy(0,2,4): x(j0,m4) - 197 y(j0,k2,m4) - x(j2,m4) >= -190
 psy(0,3,4): x(j0,m4) - 197 y(j0,k3,m4) - x(j3,m4) >= -189
 psy(0,4,4): x(j0,m4) - 197 y(j0,k4,m4) - x(j4,m4) >= -192
 psy(0,5,4): x(j0,m4) - 197 y(j0,k5,m4) - x(j5,m4) >= -193
 sum(0,4): y(j0,k1,m4) + y(j0,k2,m4) + y(j0,k3,m4) + y(j0,k4,m4)
   + y(j0,k5,m4) + y(j0,k6,m4) = 1
 psy(1,0,4): - x(j0,m4) - 197 y(j1,k0,m4) + x(j1,m4) >= -191
 psy(1,2,4): x(j1,m4) - 197 y(j1,k2,m4) - x(j2,m4) >= -190
 psy(1,3,4): x(j1,m4) - 197 y(j1,k3,m4) - x(j3,m4) >= -189
 psy(1,4,4): x(j1,m4) - 197 y(j1,k4,m4) - x(j4,m4) >= -192
 psy(1,5,4): x(j1,m4) - 197 y(j1,k5,m4) - x(j5,m4) >= -193
 sum(1,4): y(j1,k0,m4) + y(j1,k2,m4) + y(j1,k3,m4) + y(j1,k4,m4)
   + y(j1,k5,m4) + y(j1,k6,m4) = 1
 psy(2,0,4): - x(j0,m4) - 197 y(j2,k0,m4) + x(j2,m4) >= -191
 psy(2,1,4): - x(j1,m4) - 197 y(j2,k1,m4) + x(j2,m4) >= -187
 psy(2,3,4): x(j2,m4) - 197 y(j2,k3,m4) - x(j3,m4) >= -189
 psy(2,4,4): x(j2,m4) - 197 y(j2,k4,m4) - x(j4,m4) >= -192
 psy(2,5,4): x(j2,m4) - 197 y(j2,k5,m4) - x(j5,m4) >= -193
 sum(2,4): y(j2,k0,m4) + y(j2,k1,m4) + y(j2,k3,m4) + y(j2,k4,m4)
   + y(j2,k5,m4) + y(j2,k6,m4) = 1
 psy(3,0,4): - x(j0,m4) - 197 y(j3,k0,m4) + x(j3,m4) >= -191
 psy(3,1,4): - x(j1,m4) - 197 y(j3,k1,m4) + x(j3,m4) >= -187
 psy(3,2,4): - x(j2,m4) - 197 y(j3,k2,m4) + x(j3,m4) >= -190
 psy(3,4,4): x(j3,m4) - 197 y(j3,k4,m4) - x(j4,m4) >= -192
 psy(3,5,4): x(j3,m4) - 197 y(j3,k5,m4) - x(j5,m4) >= -193
 sum(3,4): y(j3,k0,m4) + y(j3,k1,m4) + y(j3,k2,m4) + y(j3,k4,m4)
   + y(j3,k5,m4) + y(j3,k6,m4) = 1
 psy(4,0,4): - x(j0,m4) - 197 y(j4,k0,m4) + x(j4,m4) >= -191
 psy(4,1,4): - x(j1,m4) - 197 y(j4,k1,m4) + x(j4,m4) >= -187
 psy(4,2,4): - x(j2,m4) - 197 y(j4,k2,m4) + x(j4,m4) >= -190
 psy(4,3,4): - x(j3,m4) - 197 y(j4,k3,m4) + x(j4,m4) >= -189
 psy(4,5,4): x(j4,m4) - 197 y(j4,k5,m4) - x(j5,m4) >= -193
 sum(4,4): y(j4,k0,m4) + y(j4,k1,m4) + y(j4,k2,m4) + y(j4,k3,m4)
   + y(j4,k5,m4) + y(j4,k6,m4) = 1
 psy(5,0,4): - x(j0,m4) - 197 y(j5,k0,m4) + x(j5,m4) >= -191
 psy(5,1,4): - x(j1,m4) - 197 y(j5,k1,m4) + x(j5,m4) >= -187
 psy(5,2,4): - x(j2,m4) - 197 y(j5,k2,m4) + x(j5,m4) >= -190
 psy(5,3,4): - x(j3,m4) - 197 y(j5,k3,m4) + x(j5,m4) >= -189
 psy(5,4,4): - x(j4,m4) - 197 y(j5,k4,m4) + x(j5,m4) >= -192
 sum(5,4): y(j5,k0,m4) + y(j5,k1,m4) + y(j5,k2,m4) + y(j5,k3,m4)
   + y(j5,k4,m4) + y(j5,k6,m4) = 1
 sum2(0,4): y(j1,k0,m4) + y(j2,k0,m4) + y(j3,k0,m4) + y(j4,k0,m4)
   + y(j5,k0,m4) <= 1
 sum2(1,4): y(j0,k1,m4) + y(j2,k1,m4) + y(j3,k1,m4) + y(j4,k1,m4)
   + y(j5,k1,m4) <= 1
 sum2(2,4): y(j0,k2,m4) + y(j1,k2,m4) + y(j3,k2,m4) + y(j4,k2,m4)
   + y(j5,k2,m4) <= 1
 sum2(3,4): y(j0,k3,m4) + y(j1,k3,m4) + y(j2,k3,m4) + y(j4,k3,m4)
   + y(j5,k3,m4) <= 1
 sum2(4,4): y(j0,k4,m4) + y(j1,k4,m4) + y(j2,k4,m4) + y(j3,k4,m4)
   + y(j5,k4,m4) <= 1
 sum2(5,4): y(j0,k5,m4) + y(j1,k5,m4) + y(j2,k5,m4) + y(j3,k5,m4)
   + y(j4,k5,m4) <= 1
 sum2(6,4): y(j0,k6,m4) + y(j1,k6,m4) + y(j2,k6,m4) + y(j3,k6,m4)
   + y(j4,k6,m4) + y(j5,k6,m4) <= 1
 psy(0,1,5): x(j0,m5) - 197 y(j0,k1,m5) - x(j1,m5) >= -187
 psy(0,2,5): x(j0,m5) - 197 y(j0,k2,m5) - x(j2,m5) >= -189
 psy(0,3,5): x(j0,m5) - 197 y(j0,k3,m5) - x(j3,m5) >= -188
 psy(0,4,5): x(j0,m5) - 197 y(j0,k4,m5) - x(j4,m5) >= -193
 psy(0,5,5): x(j0,m5) - 197 y(j0,k5,m5) - x(j5,m5) >= -188
 sum(0,5): y(j0,k1,m5) + y(j0,k2,m5) + y(j0,k3,m5) + y(j0,k4,m5)
   + y(j0,k5,m5) + y(j0,k6,m5) = 1
 psy(1,0,5): - x(j0,m5) - 197 y(j1,k0,m5) + x(j1,m5) >= -194
 psy(1,2,5): x(j1,m5) - 197 y(j1,k2,m5) - x(j2,m5) >= -189
 psy(1,3,5): x(j1,m5) - 197 y(j1,k3,m5) - x(j3,m5) >= -188
 psy(1,4,5): x(j1,m5) - 197 y(j1,k4,m5) - x(j4,m5) >= -193
 psy(1,5,5): x(j1,m5) - 197 y(j1,k5,m5) - x(j5,m5) >= -188
 sum(1,5): y(j1,k0,m5) + y(j1,k2,m5) + y(j1,k3,m5) + y(j1,k4,m5)
   + y(j1,k5,m5) + y(j1,k6,m5) = 1
 psy(2,0,5): - x(j0,m5) - 197 y(j2,k0,m5) + x(j2,m5) >= -194
 psy(2,1,5): - x(j1,m5) - 197 y(j2,k1,m5) + x(j2,m5) >= -187
 psy(2,3,5): x(j2,m5) - 197 y(j2,k3,m5) - x(j3,m5) >= -188
 psy(2,4,5): x(j2,m5) - 197 y(j2,k4,m5) - x(j4,m5) >= -193
 psy(2,5,5): x(j2,m5) - 197 y(j2,k5,m5) - x(j5,m5) >= -188
 sum(2,5): y(j2,k0,m5) + y(j2,k1,m5) + y(j2,k3,m5) + y(j2,k4,m5)
   + y(j2,k5,m5) + y(j2,k6,m5) = 1
 psy(3,0,5): - x(j0,m5) - 197 y(j3,k0,m5) + x(j3,m5) >= -194
 psy(3,1,5): - x(j1,m5) - 197 y(j3,k1,m5) + x(j3,m5) >= -187
 psy(3,2,5): - x(j2,m5) - 197 y(j3,k2,m5) + x(j3,m5) >= -189
 psy(3,4,5): x(j3,m5) - 197 y(j3,k4,m5) - x(j4,m5) >= -193
 psy(3,5,5): x(j3,m5) - 197 y(j3,k5,m5) - x(j5,m5) >= -188
 sum(3,5): y(j3,k0,m5) + y(j3,k1,m5) + y(j3,k2,m5) + y(j3,k4,m5)
   + y(j3,k5,m5) + y(j3,k6,m5) = 1
 psy(4,0,5): - x(j0,m5) - 197 y(j4,k0,m5) + x(j4,m5) >= -194
 psy(4,1,5): - x(j1,m5) - 197 y(j4,k1,m5) + x(j4,m5) >= -187
 psy(4,2,5): - x(j2,m5) - 197 y(j4,k2,m5) + x(j4,m5) >= -189
 psy(4,3,5): - x(j3,m5) - 197 y(j4,k3,m5) + x(j4,m5) >= -188
 psy(4,5,5): x(j4,m5) - 197 y(j4,k5,m5) - x(j5,m5) >= -188
 sum(4,5): y(j4,k0,m5) + y(j4,k1,m5) + y(j4,k2,m5) + y(j4,k3,m5)
   + y(j4,k5,m5) + y(j4,k6,m5) = 1
 psy(5,0,5): - x(j0,m5) - 197 y(j5,k0,m5) + x(j5,m5) >= -194
 psy(5,1,5): - x(j1,m5) - 197 y(j5,k1,m5) + x(j5,m5) >= -187
 psy(5,2,5): - x(j2,m5) - 197 y(j5,k2,m5) + x(j5,m5) >= -189
 psy(5,3,5): - x(j3,m5) - 197 y(j5,k3,m5) + x(j5,m5) >= -188
 psy(5,4,5): - x(j4,m5) - 197 y(j5,k4,m5) + x(j5,m5) >= -193
 sum(5,5): y(j5,k0,m5) + y(j5,k1,m5) + y(j5,k2,m5) + y(j5,k3,m5)
   + y(j5,k4,m5) + y(j5,k6,m5) = 1
 sum2(0,5): y(j1,k0,m5) + y(j2,k0,m5) + y(j3,k0,m5) + y(j4,k0,m5)
   + y(j5,k0,m5) <= 1
 sum2(1,5): y(j0,k1,m5) + y(j2,k1,m5) + y(j3,k1,m5) + y(j4,k1,m5)
   + y(j5,k1,m5) <= 1
 sum2(2,5): y(j0,k2,m5) + y(j1,k2,m5) + y(j3,k2,m5) + y(j4,k2,m5)
   + y(j5,k2,m5) <= 1
 sum2(3,5): y(j0,k3,m5) + y(j1,k3,m5) + y(j2,k3,m5) + y(j4,k3,m5)
   + y(j5,k3,m5) <= 1
 sum2(4,5): y(j0,k4,m5) + y(j1,k4,m5) + y(j2,k4,m5) + y(j3,k4,m5)
   + y(j5,k4,m5) <= 1
 sum2(5,5): y(j0,k5,m5) + y(j1,k5,m5) + y(j2,k5,m5) + y(j3,k5,m5)
   + y(j4,k5,m5) <= 1
 sum2(6,5): y(j0,k6,m5) + y(j1,k6,m5) + y(j2,k6,m5) + y(j3,k6,m5)
   + y(j4,k6,m5) + y(j5,k6,m5) <= 1
 makespan(0): C - x(j0,m4) >= 6
 makespan(1): C - x(j1,m3) >= 4
 makespan(2): C - x(j2,m4) >= 7
 makespan(3): C - x(j3,m5) >= 9
 makespan(4): C - x(j4,m3) >= 1
 makespan(5): C - x(j5,m2) >= 1
Bounds
Binaries
 y(j0,k1,m0) y(j1,k0,m0) y(j0,k2,m0) y(j2,k0,m0) y(j0,k3,m0) y(j3,k0,m0)
 y(j0,k4,m0) y(j4,k0,m0) y(j0,k5,m0) y(j5,k0,m0) y(j0,k6,m0) y(j1,k2,m0)
 y(j2,k1,m0) y(j1,k3,m0) y(j3,k1,m0) y(j1,k4,m0) y(j4,k1,m0) y(j1,k5,m0)
 y(j5,k1,m0) y(j1,k6,m0) y(j2,k3,m0) y(j3,k2,m0) y(j2,k4,m0) y(j4,k2,m0)
 y(j2,k5,m0) y(j5,k2,m0) y(j2,k6,m0) y(j3,k4,m0) y(j4,k3,m0) y(j3,k5,m0)
 y(j5,k3,m0) y(j3,k6,m0) y(j4,k5,m0) y(j5,k4,m0) y(j4,k6,m0) y(j5,k6,m0)
 y(j0,k1,m1) y(j1,k0,m1) y(j0,k2,m1) y(j2,k0,m1) y(j0,k3,m1) y(j3,k0,m1)
 y(j0,k4,m1) y(j4,k0,m1) y(j0,k5,m1) y(j5,k0,m1) y(j0,k6,m1) y(j1,k2,m1)
 y(j2,k1,m1) y(j1,k3,m1) y(j3,k1,m1) y(j1,k4,m1) y(j4,k1,m1) y(j1,k5,m1)
 y(j5,k1,m1) y(j1,k6,m1) y(j2,k3,m1) y(j3,k2,m1) y(j2,k4,m1) y(j4,k2,m1)
 y(j2,k5,m1) y(j5,k2,m1) y(j2,k6,m1) y(j3,k4,m1) y(j4,k3,m1) y(j3,k5,m1)
 y(j5,k3,m1) y(j3,k6,m1) y(j4,k5,m1) y(j5,k4,m1) y(j4,k6,m1) y(j5,k6,m1)
 y(j0,k1,m2) y(j1,k0,m2) y(j0,k2,m2) y(j2,k0,m2) y(j0,k3,m2) y(j3,k0,m2)
 y(j0,k4,m2) y(j4,k0,m2) y(j0,k5,m2) y(j5,k0,m2) y(j0,k6,m2) y(j1,k2,m2)
 y(j2,k1,m2) y(j1,k3,m2) y(j3,k1,m2) y(j1,k4,m2) y(j4,k1,m2) y(j1,k5,m2)
 y(j5,k1,m2) y(j1,k6,m2) y(j2,k3,m2) y(j3,k2,m2) y(j2,k4,m2) y(j4,k2,m2)
 y(j2,k5,m2) y(j5,k2,m2) y(j2,k6,m2) y(j3,k4,m2) y(j4,k3,m2) y(j3,k5,m2)
 y(j5,k3,m2) y(j3,k6,m2) y(j4,k5,m2) y(j5,k4,m2) y(j4,k6,m2) y(j5,k6,m2)
 y(j0,k1,m3) y(j1,k0,m3) y(j0,k2,m3) y(j2,k0,m3) y(j0,k3,m3) y(j3,k0,m3)
 y(j0,k4,m3) y(j4,k0,m3) y(j0,k5,m3) y(j5,k0,m3) y(j0,k6,m3) y(j1,k2,m3)
 y(j2,k1,m3) y(j1,k3,m3) y(j3,k1,m3) y(j1,k4,m3) y(j4,k1,m3) y(j1,k5,m3)
 y(j5,k1,m3) y(j1,k6,m3) y(j2,k3,m3) y(j3,k2,m3) y(j2,k4,m3) y(j4,k2,m3)
 y(j2,k5,m3) y(j5,k2,m3) y(j2,k6,m3) y(j3,k4,m3) y(j4,k3,m3) y(j3,k5,m3)
 y(j5,k3,m3) y(j3,k6,m3) y(j4,k5,m3) y(j5,k4,m3) y(j4,k6,m3) y(j5,k6,m3)
 y(j0,k1,m4) y(j1,k0,m4) y(j0,k2,m4) y(j2,k0,m4) y(j0,k3,m4) y(j3,k0,m4)
 y(j0,k4,m4) y(j4,k0,m4) y(j0,k5,m4) y(j5,k0,m4) y(j0,k6,m4) y(j1,k2,m4)
 y(j2,k1,m4) y(j1,k3,m4) y(j3,k1,m4) y(j1,k4,m4) y(j4,k1,m4) y(j1,k5,m4)
 y(j5,k1,m4) y(j1,k6,m4) y(j2,k3,m4) y(j3,k2,m4) y(j2,k4,m4) y(j4,k2,m4)
 y(j2,k5,m4) y(j5,k2,m4) y(j2,k6,m4) y(j3,k4,m4) y(j4,k3,m4) y(j3,k5,m4)
 y(j5,k3,m4) y(j3,k6,m4) y(j4,k5,m4) y(j5,k4,m4) y(j4,k6,m4) y(j5,k6,m4)
 y(j0,k1,m5) y(j1,k0,m5) y(j0,k2,m5) y(j2,k0,m5) y(j0,k3,m5) y(j3,k0,m5)
 y(j0,k4,m5) y(j4,k0,m5) y(j0,k5,m5) y(j5,k0,m5) y(j0,k6,m5) y(j1,k2,m5)
 y(j2,k1,m5) y(j1,k3,m5) y(j3,k1,m5) y(j1,k4,m5) y(j4,k1,m5) y(j1,k5,m5)
 y(j5,k1,m5) y(j1,k6,m5) y(j2,k3,m5) y(j3,k2,m5) y(j2,k4,m5) y(j4,k2,m5)
 y(j2,k5,m5) y(j5,k2,m5) y(j2,k6,m5) y(j3,k4,m5) y(j4,k3,m5) y(j3,k5,m5)
 y(j5,k3,m5) y(j3,k6,m5) y(j4,k5,m5) y(j5,k4,m5) y(j4,k6,m5) y(j5,k6,m5)
Generals
 C x(j0,m0) x(j1,m0) x(j2,m0) x(j3,m0) x(j4,m0) x(j5,m0) x(j0,m1) x(j1,m1)
 x(j2,m1) x(j3,m1) x(j4,m1) x(j5,m1) x(j0,m2) x(j1,m2) x(j2,m2) x(j3,m2)
 x(j4,m2) x(j5,m2) x(j0,m3) x(j1,m3) x(j2,m3) x(j3,m3) x(j4,m3) x(j5,m3)
 x(j0,m4) x(j1,m4) x(j2,m4) x(j3,m4) x(j4,m4) x(j5,m4) x(j0,m5) x(j1,m5)
 x(j2,m5) x(j3,m5) x(j4,m5) x(j5,m5)
End
