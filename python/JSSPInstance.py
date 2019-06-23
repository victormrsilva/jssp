#JSSPInstance.py
class JSSPInstance:
    def __init__(self, instanceName):
        self.instanceName = instanceName
        with open(instanceName,"r") as instance:
            line = instance.readline()
            self.n, self.m = (int(x) for x in line.split())
            self.machines = [[0]*self.m for i in range(self.n)]
            self.times = [[0]*self.m for i in range(self.n)]
            self.e = [[0]*self.m for i in range(self.n)]
            self.f = [[0]*self.m for i in range(self.n)]
            self.K = 0
            
            #read instance file
            for j in range(self.n):
                line = instance.readline()
                line = line.split()
                i = 0
                while i < self.m:
                    value = int(line[2*i])
                    self.machines[j][i] = value
                    value = int(line[2*i+1])
                    self.times[j][int(line[2*i])] = value
                    self.K += value
                    
                    i = i+1

            #fill e and f matrix
            for j in range(self.n):
                for i in range(self.m):
                    for aux in range(0,i):
                        self.e[j][self.machines[j][i]] += self.times[j][self.machines[j][aux]]

                    for aux in range(i+1,self.m):
                        self.f[j][self.machines[j][i]] += self.times[j][self.machines[j][aux]]

    
    def print(self):
        print("Machines: ", self.m, "Jobs: ", self.n)
        print("Order: ")
        for j in range(self.n):
            print("Job ",j+1,": ", end='')
            for i in range(self.m):
                print("machine",self.machines[j][i]+1, "(", self.times[j][self.machines[j][i]], ") ", end='')
            print()
        print("K: ",self.K)
        for j in range(self.n):
            print("Job ",j+1,", e: ", end='')
            for i in range(self.m):
                print("machine ",i+1, "(", self.e[j][i], ") ", end='')
            print()
        for j in range(self.n):
            print("Job ",j+1,", f: ", end='')
            for i in range(self.m):
                print("machine ",i+1, "(", self.f[j][i], ") ", end='')
            print()


