#flow.py
from mip.model import Model
import sys
from time import process_time
from instance import Instance

class Flow:
  def __init__(self, instance):
    self.instance = instance

  def constructProblem(self):
    self.instance.print()
    self.model = Model('flow')
    self.fIdx = [[self.model.add_var('f({},{})'.format(i+1,t), var_type='B') for i in range(self.instance.m())] for t in range(self.instance.h())]
    self.xIdx = [[[self.model.add_var('x({},{},{})'.format(i+1,j+1,t)) for t in range(self.instance.est(i,j),self.instance.lst(i,j)+1)] for i in range(self.instance.m())] for j in range(self.instance.n())]
    self.eIdx = [[[self.model.add_var('e({},{},{})'.format(i+1,j+1,t)) for t in range(self.instance.est(i,j),self.instance.lst(i,j)+1)] for i in range(self.instance.m())] for j in range(self.instance.n())]
    self.cIdx = self.model.add_var('C', var_type='I')
    
    for aux in self.fIdx:
      for f in aux:
        print(f.name," ",end='')
      print()
    for aux in self.xIdx:
      for aux2 in aux:
        for x in aux2:
          print(x.name," ",end='')
        print()
    for aux in self.xIdx:
      for aux2 in aux:
        for e in aux2:
          print(e.name," ",end='')
        print()




