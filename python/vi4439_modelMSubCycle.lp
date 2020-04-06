\ Model compact
\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  C
Subject To
 ord(0,1): - x(j0,m0) + x(j0,m3) >= 10
 ord(0,2): x(j0,m2) - x(j0,m3) >= 1
 ord(0,3): x(j0,m1) - x(j0,m2) >= 6
 ord(1,1): - x(j1,m0) + x(j1,m3) >= 1
 ord(1,2): x(j1,m1) - x(j1,m3) >= 7
 ord(1,3): - x(j1,m1) + x(j1,m2) >= 5
 ord(2,1): x(j2,m2) - x(j2,m3) >= 9
 ord(2,2): x(j2,m1) - x(j2,m2) >= 3
 ord(2,3): x(j2,m0) - x(j2,m1) >= 7
 ord(3,1): - x(j3,m1) + x(j3,m3) >= 7
 ord(3,2): x(j3,m0) - x(j3,m3) >= 3
 ord(3,3): - x(j3,m0) + x(j3,m2) >= 9
 psy(0,1,0): x(j0,m0) - 94 y(j0,k1,m0) - x(j1,m0) >= -93
 psy(0,2,0): x(j0,m0) - 94 y(j0,k2,m0) - x(j2,m0) >= -93
 psy(0,3,0): x(j0,m0) - 94 y(j0,k3,m0) - x(j3,m0) >= -85
 sum(0,0): y(j0,k1,m0) + y(j0,k2,m0) + y(j0,k3,m0) + y(j0,k4,m0) = 1
 psy(1,0,0): - x(j0,m0) - 94 y(j1,k0,m0) + x(j1,m0) >= -84
 psy(1,2,0): x(j1,m0) - 94 y(j1,k2,m0) - x(j2,m0) >= -93
 psy(1,3,0): x(j1,m0) - 94 y(j1,k3,m0) - x(j3,m0) >= -85
 sum(1,0): y(j1,k0,m0) + y(j1,k2,m0) + y(j1,k3,m0) + y(j1,k4,m0) = 1
 psy(2,0,0): - x(j0,m0) - 94 y(j2,k0,m0) + x(j2,m0) >= -84
 psy(2,1,0): - x(j1,m0) - 94 y(j2,k1,m0) + x(j2,m0) >= -93
 psy(2,3,0): x(j2,m0) - 94 y(j2,k3,m0) - x(j3,m0) >= -85
 sum(2,0): y(j2,k0,m0) + y(j2,k1,m0) + y(j2,k3,m0) + y(j2,k4,m0) = 1
 psy(3,0,0): - x(j0,m0) - 94 y(j3,k0,m0) + x(j3,m0) >= -84
 psy(3,1,0): - x(j1,m0) - 94 y(j3,k1,m0) + x(j3,m0) >= -93
 psy(3,2,0): - x(j2,m0) - 94 y(j3,k2,m0) + x(j3,m0) >= -93
 sum(3,0): y(j3,k0,m0) + y(j3,k1,m0) + y(j3,k2,m0) + y(j3,k4,m0) = 1
 sum2(0,0): y(j1,k0,m0) + y(j2,k0,m0) + y(j3,k0,m0) <= 1
 sum2(1,0): y(j0,k1,m0) + y(j2,k1,m0) + y(j3,k1,m0) <= 1
 sum2(2,0): y(j0,k2,m0) + y(j1,k2,m0) + y(j3,k2,m0) <= 1
 sum2(3,0): y(j0,k3,m0) + y(j1,k3,m0) + y(j2,k3,m0) <= 1
 sum2(4,0): y(j0,k4,m0) + y(j1,k4,m0) + y(j2,k4,m0) + y(j3,k4,m0) <= 1
 psy(0,1,1): x(j0,m1) - 94 y(j0,k1,m1) - x(j1,m1) >= -89
 psy(0,2,1): x(j0,m1) - 94 y(j0,k2,m1) - x(j2,m1) >= -87
 psy(0,3,1): x(j0,m1) - 94 y(j0,k3,m1) - x(j3,m1) >= -87
 sum(0,1): y(j0,k1,m1) + y(j0,k2,m1) + y(j0,k3,m1) + y(j0,k4,m1) = 1
 psy(1,0,1): - x(j0,m1) - 94 y(j1,k0,m1) + x(j1,m1) >= -87
 psy(1,2,1): x(j1,m1) - 94 y(j1,k2,m1) - x(j2,m1) >= -87
 psy(1,3,1): x(j1,m1) - 94 y(j1,k3,m1) - x(j3,m1) >= -87
 sum(1,1): y(j1,k0,m1) + y(j1,k2,m1) + y(j1,k3,m1) + y(j1,k4,m1) = 1
 psy(2,0,1): - x(j0,m1) - 94 y(j2,k0,m1) + x(j2,m1) >= -87
 psy(2,1,1): - x(j1,m1) - 94 y(j2,k1,m1) + x(j2,m1) >= -89
 psy(2,3,1): x(j2,m1) - 94 y(j2,k3,m1) - x(j3,m1) >= -87
 sum(2,1): y(j2,k0,m1) + y(j2,k1,m1) + y(j2,k3,m1) + y(j2,k4,m1) = 1
 psy(3,0,1): - x(j0,m1) - 94 y(j3,k0,m1) + x(j3,m1) >= -87
 psy(3,1,1): - x(j1,m1) - 94 y(j3,k1,m1) + x(j3,m1) >= -89
 psy(3,2,1): - x(j2,m1) - 94 y(j3,k2,m1) + x(j3,m1) >= -87
 sum(3,1): y(j3,k0,m1) + y(j3,k1,m1) + y(j3,k2,m1) + y(j3,k4,m1) = 1
 sum2(0,1): y(j1,k0,m1) + y(j2,k0,m1) + y(j3,k0,m1) <= 1
 sum2(1,1): y(j0,k1,m1) + y(j2,k1,m1) + y(j3,k1,m1) <= 1
 sum2(2,1): y(j0,k2,m1) + y(j1,k2,m1) + y(j3,k2,m1) <= 1
 sum2(3,1): y(j0,k3,m1) + y(j1,k3,m1) + y(j2,k3,m1) <= 1
 sum2(4,1): y(j0,k4,m1) + y(j1,k4,m1) + y(j2,k4,m1) + y(j3,k4,m1) <= 1
 psy(0,1,2): x(j0,m2) - 94 y(j0,k1,m2) - x(j1,m2) >= -86
 psy(0,2,2): x(j0,m2) - 94 y(j0,k2,m2) - x(j2,m2) >= -91
 psy(0,3,2): x(j0,m2) - 94 y(j0,k3,m2) - x(j3,m2) >= -84
 sum(0,2): y(j0,k1,m2) + y(j0,k2,m2) + y(j0,k3,m2) + y(j0,k4,m2) = 1
 psy(1,0,2): - x(j0,m2) - 94 y(j1,k0,m2) + x(j1,m2) >= -88
 psy(1,2,2): x(j1,m2) - 94 y(j1,k2,m2) - x(j2,m2) >= -91
 psy(1,3,2): x(j1,m2) - 94 y(j1,k3,m2) - x(j3,m2) >= -84
 sum(1,2): y(j1,k0,m2) + y(j1,k2,m2) + y(j1,k3,m2) + y(j1,k4,m2) = 1
 psy(2,0,2): - x(j0,m2) - 94 y(j2,k0,m2) + x(j2,m2) >= -88
 psy(2,1,2): - x(j1,m2) - 94 y(j2,k1,m2) + x(j2,m2) >= -86
 psy(2,3,2): x(j2,m2) - 94 y(j2,k3,m2) - x(j3,m2) >= -84
 sum(2,2): y(j2,k0,m2) + y(j2,k1,m2) + y(j2,k3,m2) + y(j2,k4,m2) = 1
 psy(3,0,2): - x(j0,m2) - 94 y(j3,k0,m2) + x(j3,m2) >= -88
 psy(3,1,2): - x(j1,m2) - 94 y(j3,k1,m2) + x(j3,m2) >= -86
 psy(3,2,2): - x(j2,m2) - 94 y(j3,k2,m2) + x(j3,m2) >= -91
 sum(3,2): y(j3,k0,m2) + y(j3,k1,m2) + y(j3,k2,m2) + y(j3,k4,m2) = 1
 sum2(0,2): y(j1,k0,m2) + y(j2,k0,m2) + y(j3,k0,m2) <= 1
 sum2(1,2): y(j0,k1,m2) + y(j2,k1,m2) + y(j3,k1,m2) <= 1
 sum2(2,2): y(j0,k2,m2) + y(j1,k2,m2) + y(j3,k2,m2) <= 1
 sum2(3,2): y(j0,k3,m2) + y(j1,k3,m2) + y(j2,k3,m2) <= 1
 sum2(4,2): y(j0,k4,m2) + y(j1,k4,m2) + y(j2,k4,m2) + y(j3,k4,m2) <= 1
 psy(0,1,3): x(j0,m3) - 94 y(j0,k1,m3) - x(j1,m3) >= -87
 psy(0,2,3): x(j0,m3) - 94 y(j0,k2,m3) - x(j2,m3) >= -85
 psy(0,3,3): x(j0,m3) - 94 y(j0,k3,m3) - x(j3,m3) >= -91
 sum(0,3): y(j0,k1,m3) + y(j0,k2,m3) + y(j0,k3,m3) + y(j0,k4,m3) = 1
 psy(1,0,3): - x(j0,m3) - 94 y(j1,k0,m3) + x(j1,m3) >= -93
 psy(1,2,3): x(j1,m3) - 94 y(j1,k2,m3) - x(j2,m3) >= -85
 psy(1,3,3): x(j1,m3) - 94 y(j1,k3,m3) - x(j3,m3) >= -91
 sum(1,3): y(j1,k0,m3) + y(j1,k2,m3) + y(j1,k3,m3) + y(j1,k4,m3) = 1
 psy(2,0,3): - x(j0,m3) - 94 y(j2,k0,m3) + x(j2,m3) >= -93
 psy(2,1,3): - x(j1,m3) - 94 y(j2,k1,m3) + x(j2,m3) >= -87
 psy(2,3,3): x(j2,m3) - 94 y(j2,k3,m3) - x(j3,m3) >= -91
 sum(2,3): y(j2,k0,m3) + y(j2,k1,m3) + y(j2,k3,m3) + y(j2,k4,m3) = 1
 psy(3,0,3): - x(j0,m3) - 94 y(j3,k0,m3) + x(j3,m3) >= -93
 psy(3,1,3): - x(j1,m3) - 94 y(j3,k1,m3) + x(j3,m3) >= -87
 psy(3,2,3): - x(j2,m3) - 94 y(j3,k2,m3) + x(j3,m3) >= -85
 sum(3,3): y(j3,k0,m3) + y(j3,k1,m3) + y(j3,k2,m3) + y(j3,k4,m3) = 1
 sum2(0,3): y(j1,k0,m3) + y(j2,k0,m3) + y(j3,k0,m3) <= 1
 sum2(1,3): y(j0,k1,m3) + y(j2,k1,m3) + y(j3,k1,m3) <= 1
 sum2(2,3): y(j0,k2,m3) + y(j1,k2,m3) + y(j3,k2,m3) <= 1
 sum2(3,3): y(j0,k3,m3) + y(j1,k3,m3) + y(j2,k3,m3) <= 1
 sum2(4,3): y(j0,k4,m3) + y(j1,k4,m3) + y(j2,k4,m3) + y(j3,k4,m3) <= 1
 makespan(0): C - x(j0,m1) >= 7
 makespan(1): C - x(j1,m2) >= 8
 makespan(2): C - x(j2,m0) >= 1
 makespan(3): C - x(j3,m2) >= 10
Bounds
Binaries
 y(j0,k1,m0) y(j1,k0,m0) y(j0,k2,m0) y(j2,k0,m0) y(j0,k3,m0) y(j3,k0,m0)
 y(j0,k4,m0) y(j1,k2,m0) y(j2,k1,m0) y(j1,k3,m0) y(j3,k1,m0) y(j1,k4,m0)
 y(j2,k3,m0) y(j3,k2,m0) y(j2,k4,m0) y(j3,k4,m0) y(j0,k1,m1) y(j1,k0,m1)
 y(j0,k2,m1) y(j2,k0,m1) y(j0,k3,m1) y(j3,k0,m1) y(j0,k4,m1) y(j1,k2,m1)
 y(j2,k1,m1) y(j1,k3,m1) y(j3,k1,m1) y(j1,k4,m1) y(j2,k3,m1) y(j3,k2,m1)
 y(j2,k4,m1) y(j3,k4,m1) y(j0,k1,m2) y(j1,k0,m2) y(j0,k2,m2) y(j2,k0,m2)
 y(j0,k3,m2) y(j3,k0,m2) y(j0,k4,m2) y(j1,k2,m2) y(j2,k1,m2) y(j1,k3,m2)
 y(j3,k1,m2) y(j1,k4,m2) y(j2,k3,m2) y(j3,k2,m2) y(j2,k4,m2) y(j3,k4,m2)
 y(j0,k1,m3) y(j1,k0,m3) y(j0,k2,m3) y(j2,k0,m3) y(j0,k3,m3) y(j3,k0,m3)
 y(j0,k4,m3) y(j1,k2,m3) y(j2,k1,m3) y(j1,k3,m3) y(j3,k1,m3) y(j1,k4,m3)
 y(j2,k3,m3) y(j3,k2,m3) y(j2,k4,m3) y(j3,k4,m3)
Generals
 C x(j0,m0) x(j1,m0) x(j2,m0) x(j3,m0) x(j0,m1) x(j1,m1) x(j2,m1) x(j3,m1)
 x(j0,m2) x(j1,m2) x(j2,m2) x(j3,m2) x(j0,m3) x(j1,m3) x(j2,m3) x(j3,m3)
End
