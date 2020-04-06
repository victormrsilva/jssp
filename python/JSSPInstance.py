#JSSPInstance.py
class JSSPInstance:
    def __init__(self, instancename, H):
        self.instancename = instancename
        with open(instancename, "r") as instance:
            line = instance.readline()
            self.n, self.m = (int(x) for x in line.split())
            self.machines = [[0]*self.m for i in range(self.n)]
            self.o = [[0]*self.m for i in range(self.n)]
            self.times = [[0]*self.m for i in range(self.n)]
            self.e = [[0]*self.m for i in range(self.n)]
            self.f = [[0]*self.m for i in range(self.n)]
            if H > 0:
                self.K = H
            else:
                self.K = 0

            # read instance file
            for j in range(self.n):
                line = instance.readline()
                line = line.split()
                i = 0
                while i < self.m:
                    value = int(line[2*i])
                    self.machines[j][i] = value
                    self.o[j][value] = i
                    value = int(line[2*i+1])
                    self.times[j][int(line[2*i])] = value
                    if H < 0:
                        self.K += value

                    i = i+1

            # fill e and f matrix
            for j in range(self.n):
                for i in range(self.m):
                    for aux in range(0,i):
                        self.e[j][self.machines[j][i]] += self.times[j][self.machines[j][aux]]

                    for aux in range(i+1, self.m):
                        self.f[j][self.machines[j][i]] += self.times[j][self.machines[j][aux]]

            self.est = [[0]*self.m for i in range(self.n)]
            self.lst = [[0]*self.m for i in range(self.n)]
            self.distances = [[[0]*self.m for i in range(self.m)] for j in range(self.n)]

            for j in range(self.n):
                t = 0
                for i in range(self.m):
                    m = self.machines[j][i]
                    self.est[j][m] = t
                    t = t + self.times[j][m]

            for j in range(self.n):
                t = self.K
                for i in range(self.m-1, -1, -1):
                    m = self.machines[j][i]
                    t = t - self.times[j][m]
                    self.lst[j][m] = t

            for j in range(self.n):
                for i in range(self.m):
                    ma = self.machines[j][i]
                    for k in range(i+1,self.m):
                        mb = self.machines[j][k]
                        dur = self.times[j][ma]
                        for aux in range(i+1,k):
                            dur = dur + self.times[j][self.machines[j][aux]]
                        self.distances[j][ma][mb] += dur

    def print(self):
        print("Machines: ", self.m, "Jobs: ", self.n)
        print("Order: ")
        for j in range(self.n):
            print("Job ", j, ": ", end='')
            for i in range(self.m):
                print("machine", self.machines[j][i], "(", self.times[j][self.machines[j][i]], ") ", end='')
            print()
        print("K: ", self.K)
        for j in range(self.n):
            print("Job ", j, ", e: ", end='')
            for i in range(self.m):
                print("machine ", i, "(", self.e[j][i], ") ", end='')
            print()
        for j in range(self.n):
            print('Job ', j, ', f: ', end='')
            for i in range(self.m):
                print('machine ', i, '(', self.f[j][i], ") ", end='')
            print()
        print("EST and LST:")
        for j in range(self.n):
            for i in range(self.m):
                print("Job: {} Machine: {} EST: {} LST: {}".format(j, i, self.est[j][i],self.lst[j][i]))
        print("Order:")
        for j in range(self.n):
            for i in range(self.m):
                print("Job: {} Machine: {} = {} ".format(j, i, self.o[j][i]))
        print("Distances:")
        for j in range(self.n):
            for i in range(self.m):
                for k in range(i+1,self.m):
                  print("Job {}, Machine {} -> Machine {}: distance {}".format(j, self.machines[j][i], self.machines[j][k], self.distances[j][self.machines[j][i]][self.machines[j][k]]))


