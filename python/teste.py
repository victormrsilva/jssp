from mip.model import *
p = [10, 13, 18, 31,  7, 15]
w = [11, 15, 20, 35, 10, 33]
c = 40
n = len(w)
m = Model('knapsack', MAXIMIZE,solver_name="gurobi")
x = [m.add_var(var_type='B') for i in range(n)]
m += xsum(p[i]*x[i] for i in range(n) )
m += xsum(w[i]*x[i] for i in range(n) ) <= c
m.optimize()
selected=[i for i in range(n) if x[i].x>=0.99]
print('selected items: {}'.format(selected))
