from time import time
from compact import Compact
from mip import OptimizationStatus
from node import Node
import collections

class Branch:

    def branch(self, compact):
        bestUB = compact.instance.K

        # adding optimize relaxed problem and adding cuts
        node = Node(compact.model, compact.instance)
        # compact.relax()
        # compact.model.write('teste1.lp')

        node.optNode0()
        node.depth = 0
        # node.mip.write('testenode.lp')
        # input('testenode')

        L = collections.deque([node])

        while len(L) != 0:

            p = L.pop()

            lb = p.mip.objective_value
            print('progress: LB = {} - bestUB = {} - gap = {}%'.format(lb, bestUB,
                                                                       round((bestUB - lb) / bestUB * 100, 2)))

            dict = {}
            # input('iniciando dicionario')
            for i in range(p.instance.m):
                for j in range(p.instance.n):
                    for k in range(j + 1, p.instance.n):
                        s = (i, j, k)
                        dist = min(abs(p.y[i][j][k].x - 0.5), abs(p.y[i][k][j].x - 0.5))
                        valuejk = p.valueOrder(j, k, i)
                        valuekj = p.valueOrder(k, j, i)
                        if max(valuejk, valuekj) == 0:
                            print('max(valuejk, valuekj) == 0 ; j = {} k = {} i = {} valuejk = {} valuekj = {}'.format(j, k, i, valuejk, valuekj))
                            continue
                        value = abs(valuejk - valuekj) / (2*max(valuejk, valuekj))  # ficar valor entre 0 e 0,5
                        # input('{} {} {} {} {} {} {}'.format(j, k, i, valuejk, valuekj, value, dist))
                        dist = dist + value
                        if value < 0.75:
                            dict[s] = dist
            # print(dict)
            a = sorted(dict.items(), key=lambda l: l[1])
            # print('ordenado')
            # print(a)
            eps = 0.00000001
            # input(len(a))
            qtd = 0
            start_iter = time()
            for key, value in a:
                i = key[0]
                j = key[1]
                k = key[2]
                qtd = qtd + 1
                p1 = Node(p.mip, p.instance)  # j before k in machine i
                p1.depth = p.depth + 1
                p2 = Node(p.mip, p.instance)  # k before j in machine i
                p2.depth = p.depth + 1
                # input('atualizando nós')
                start = time()
                p1.updateNode(j, k, i)
                end = time()
                print('Time elapsed relaxation p1: {}s ; '.format(round(end - start, 2)), end='')
                # input('teste update node')
                # p1.printSolution()
                
                p1.mip.write('testep1depth{}.lp'.format(p1.depth))
                if p1.mip.status == OptimizationStatus.OPTIMAL or p1.mip.status == OptimizationStatus.FEASIBLE:  # optimal or feasible solution found
                    print('obj p1: {}'.format(p1.mip.objective_value), end='')
                    if self.is_whole(p1.mip.objective_value, eps):
                        print(' obj integer', end='')
                        if self.integerSol(p1):
                            # verificar se é melhor e então atualizar
                            print(' ; solução inteira encontrada: {}'.format(p1.mip.objective_value))
                            bestUB = p1.mip.objective_value
                        else:
                            if p1.mip.objective_value - (bestUB - 1 + eps) <= eps:
                                L.append(p1)
                    else:
                        print(' obj float')
                        if p1.mip.objective_value - (bestUB - 1 + eps) <= eps:
                            L.append(p1)
                else:
                    print(' infesiable')
                
                start = time()
                p2.updateNode(k, j, i)
                end = time()
                print('Time elapsed relaxation p2: {}s ; '.format(round(end - start, 2)), end='')

                p2.mip.write('testep2depth{}.lp'.format(p2.depth))
                if p2.mip.status == OptimizationStatus.OPTIMAL or p2.mip.status == OptimizationStatus.FEASIBLE:  # feasible solution found
                    print('obj p2: {}'.format(p2.mip.objective_value), end='')
                    if self.is_whole(p2.mip.objective_value, eps):
                        print(' obj integer', end='')
                        if self.integerSol(p2):
                            print('solução inteira encontrada: {}'.format(p2.mip.objective_value))
                            bestUB = p2.mip.objective_value
                        else:
                            if p2.mip.objective_value - (bestUB - 1 + eps) <= eps:
                                L.append(p2)
                    else:
                        print(' obj float')
                        if p2.mip.objective_value - (bestUB - 1 + eps) <= eps:
                            L.append(p2)

                # remover de L os que estão violados (bound >= bestUB)
                for index in range(len(L) - 1, -1, -1):
                    if L[index].mip.objective_value - bestUB > eps:
                        del L[index]

                print('quantidade de nós vistos no depth {}: {} de {}'.format(p.depth, qtd, len(a)))
                print('list size: {}'.format(len(L)))
            end_iter = time()
            print('Time elapsed for iteration: {}s'.format(round(end_iter - start_iter, 2)))
            print('new depth')
    def integerSol(self, p):
        for i in range(p.instance.m):
            for j in range(p.instance.n):
                if not self.is_whole(p.x[j][i].x, 0.00000001):
                    print(' ; var not integer: {} value {}'.format(p.x[j][i].name, p.x[j][i].x))
                    return False
        return True

    def is_whole(self, f, eps):
        return abs(f - round(f)) < abs(eps)
