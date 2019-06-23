from sys import argv
from JSSPInstance import JSSPInstance
from compact import Compact

inst = JSSPInstance(argv[1])

compact = Compact(inst)

compact.constructProblem()

# model = Model('jssp')

# c = model.add_var(var_type=INTEGER, name="C")
# x = [[model.add_var(var_type=INTEGER, name='x({},{})'.format(j+1,i+1)) for i in range(inst.m)] for j in range(inst.n)]
# y = [[[model.add_var(var_type=BINARY, name='y({},{},{})'.format(j+1,k+1,i+1)) for i in range(inst.m)] for k in range(inst.n)] for j in range(inst.n)]

# model.objective = c

# #constraints (2)
# for j in range(inst.n):
#     for i in range(1,inst.m):
#         model += x[j][inst.machines[j][i]] - x[j][inst.machines[j][i-1]] >= inst.times[j][inst.machines[j][i-1]] , 'ord({},{})'.format(j+1,i+1)

# #constraints (3-4)
# for j in range(inst.n):
#     for k in range(inst.n):
#         if k != j:
#             for i in range(inst.m):
#                 model += x[j][i] - x[k][i] + inst.K*y[j][k][i] >= inst.times[k][i], 'phi({},{},{})'.format(j+1,k+1,i+1)
#                 model += -x[j][i] + x[k][i] - inst.K*y[j][k][i] >= inst.times[j][i] - inst.K, 'psy({},{},{})'.format(j+1,k+1,i+1)
            
# #constraints (5)
# for j in range(inst.n):
#     model += c - x[j][inst.machines[j][inst.m - 1]] >= inst.times[j][inst.machines[j][inst.m - 1]], 'makespan({})'.format(j+1)

# model.optimize()

# #printing results
# print("C: ", c.x)
# for j in range(inst.n):
#     for i in range(inst.m):
#         print('x({},{}) = {} '.format(j+1,i+1,x[j][i].x), end='')
#     print()

