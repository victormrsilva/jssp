#instance.py
class Instance:
  def __init__(self, instanceName, maxTime):
    self.instanceName = instanceName
    self.maxTime = maxTime
    with open(instanceName,"r") as instance:
      line = instance.readline()
      self.n_, self.m_ = (int(x) for x in line.split())
      self.machines = [[0]*self.m_ for i in range(self.n_)]
      self.times = [[0]*self.m_ for i in range(self.n_)]
      
      for i in range(self.n_):
        line = instance.readline()
        line = line.split()
        j = 0
        while j < self.m_:
          value = int(line[2*j])
          self.machines[i][j] = value
          value = int(line[2*j+1])
          self.times[i][int(line[2*j])] = value
          j = j+1
      self.est_ = [[0]*self.m_ for i in range(self.n_)]
      self.lst_ = [[0]*self.m_ for i in range(self.n_)]
      self.distances = [[[0]*self.m_ for i in range(self.m_)] for j in range(self.n_)]

      for j in range(self.n_):
        t = 0
        for i in range(self.m_):
          m = self.machines[j][i]
          self.est_[j][m] = t
          t = t + self.times[j][m]

      for j in range(self.n_):
        t = self.maxTime
        for i in range(self.m_-1, -1, -1):
          m = self.machines[j][i]
          t = t - self.times[j][m]
          self.lst_[j][m] = t

      for j in range(self.n_):
        for i in range(self.m_):
          m_a = self.machines[j][i]
          for k in range(i+1,self.m_):
            m_b = self.machines[j][k]
            dur = self.times[j][m_a]
            for aux in range(i+1,k):
              dur = dur + self.times[j][self.machines[j][aux]]
            self.distances[j][m_a][m_b] += dur
  
  def print(self):
    print("Instance: ",self.instanceName)
    print("Machines: ", self.m_, " Jobs: ", self.n_, " Time Horizon: ", self.maxTime)
    print("Order: ")
    for i in range(self.n_):
      print("Job ",i+1,": ", end='')
      for j in range(self.m_):
        print("machine",self.machines[i][j]+1, "(", self.times[i][self.machines[i][j]], ") ", end='')
      print()

    print("EST and LST:")
    for j in range(self.n_):
      for i in range(self.m_):
        print("Job: {} Machine: {} EST: {} LST: {}".format(j+1,i+1,self.est_[j][i],self.lst_[j][i]))
    print("Distances:")
    for j in range(self.n_):
      for i in range(self.m_):
        for k in range(i+1,self.m_):
          print("Job {}, Machine {} -> Machine {}: distance {}".format(j+1,self.machines[j][i]+1,self.machines[j][k]+1,self.distances[j][self.machines[j][i]][self.machines[j][k]]))

  def h(self):
    return self.maxTime

  def getMachine(self,job,operation):
    return self.machines[job][operation]
  
  def getTime(self,job,machine):
    return self.times[job][machine]

  def est(self,job,machine):
    return self.est_[job][machine]    

  def lst(self,job,machine):
    return self.lst_[job][machine]    

  def getInstanceName(self):
    return self.instanceName

  def m(self):
    return self.m_

  def n(self):
    return self.n_

  
    