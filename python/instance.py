#instance.py
class Instance:
  def __init__(self, instanceName, maxTime):
    self.instanceName = instanceName
    self.maxTime = maxTime
    with open(instanceName,"r") as instance:
      line = instance.readline()
      self.n, self.m = (int(x) for x in line.split())
      self.machine = [[0]*self.m for i in range(self.n)]
      self.times = [[0]*self.m for i in range(self.n)]
      
      for i in range(self.n):
        line = instance.readline()
        line = line.split()
        j = 0
        while j < self.m:
          value = int(line[2*j])
          self.machine[i][j] = value
          value = int(line[2*j+1])
          self.times[i][int(line[2*j])] = value
          j = j+1

  
  def print(self):
    print("Instance: ",self.instanceName)
    print("Machines: ", self.m, " Tasks: ", self.n, " Time Horizon: ", self.maxTime)
    print("Order: ")
    for i in range(self.n):
      print("task ",i+1)
      for j in range(self.m):
        print("machine ",self.machine[i][j]+1, " (", self.times[i][self.machine[i][j]], ") ", end='')
      print()

  def getMaxTime(self):
    return self.maxTime

  def getMachine(self,task,operation):
    return self.machine[task][operation]
  
  def getTime(self,task,machine):
    return self.times[task][machine]

  def getInstanceName(self):
    return self.instanceName

    
    