from itertools import permutations, combinations

import numpy as np
from mip.model import Model, xsum
from mip.constants import INTEGER, BINARY, CONTINUOUS


# s = list(range(0, 25))
# print(s)

# p = []
# for i in range(2, 6):
#     seq = permutations(s, i)
#     p += list(seq)


# x = []
# c = 3
# tam = 5
# s = list(range(tam))
# print(s)
# s_aux = list(range(tam))
# perms = np.asarray(list(permutations(s_aux, c)))
# perms = [np.asarray(s) for s in perms]
# p = []
# for i in range(2, c+1):
#     seq = permutations(s, i)
#     p += list(seq)
# 
# print(p)
# 
# x = [[] for i in range(c)]
# 
# x[0] = ['x({})'.format(i) for i in range(tam)]
# x[1] = [['x({},{})'.format(a, b) for a in range(tam)] for b in range(tam)]
# x[2] = [[['x({},{},{})'.format(a, b, c) for a in range(tam)] for b in range(tam)] for c in range(tam)]
# 
# for a in range(tam):
#     for b in range(a+1, tam):
#         x[1][b][a] = x[1][a][b]
#         for c in range(b+1, tam):
#             x[2][a][c][b] = x[2][a][b][c]
#             x[2][c][a][b] = x[2][a][b][c]
#             x[2][b][c][a] = x[2][a][b][c]
#             x[2][b][a][c] = x[2][a][b][c]
# 
# for a in range(tam):
#     print('({}) = {})'.format(a, x[0][a]))
#     for b in range(tam):
#         print('({},{}) = {})'.format(a, b, x[1][a][b]))
#         for c in range(tam):
#             print('({},{},{}) = {})'.format(a, b, c, x[2][a][b][c]))


# print('Best triangle cuts')

# y = [[0, 0.7, 0,2], [0.7, 0, 0.8], [0.6, 0.2, 0]]
# m = Model('')
# m.verbose = 0
# xij = [[m.add_var(var_type=BINARY, lb=0, name='xi_xj({},{})'.format(i, j)) for i in range(3)]
#        for j in range(3)]
#
# var = - xsum(xij[i][j] for i in range(3) for j in range(3) if i != j)
# print(var)
# # input()
# m.objective = var
# # print(var)
#
# for i in range(3):
#     m += xsum(xij[i][j] for j in range(3) if j != i) - xsum(xij[j][i] for j in range(3) if j != i) == 0, 'flow({})'.format(i)
#
# m += xsum(y[i][j] * xij[i][j] for i in range(3) for j in range(3) if j != i
#           ) - xsum(xij[i][j] for i in range(3) for j in range(3) if j != i) >= -1, 'triangle_cut'
#
# m.write('best_model.lp')
# # input('modelo pronto')
# m.optimize()
# print('{}'.format((round(m.objective_value, 2))))
# for i in range(3):
#     for j in range(3):
#         print('x[{}][{}] = {}'.format(i, j, round(xij[i][j].x, 2)))
#
# S = []
# for i in range(3):
#     for j in range(3):
#         if xij[i][j].x > 0.99999:
#             S.append(y[i][j])
#
# var = xsum(i for i in S) <= len(S) - 1
# print(var)
# input()




x = [0, 0, 19, 10]

print(x)

for sizeS in range(1, 4):
    comb = combinations(list(range(0, 4)), sizeS + 1)
    comb = list(comb)
    dict = {}
    for s in range(len(comb)):
        dist = 0
        S = comb[s]
        for i in list(range(0, sizeS+1)):
            for j in list(range(i+1, sizeS+1)):
                dist += abs(x[S[i]] - x[S[j]])
        print('{}: {}'.format(S, dist))
        dict[s] = dist
    print(dict)
    print(sorted(dict.items(), key=lambda l: l[1]))
    input()


exit()
