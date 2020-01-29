import string
from time import time
from math import modf
import re
from compact import Compact
from mip import OptimizationStatus
from node import Node
import collections

class Branch:

    def branch(self, compact):
        print('time\t;node\t;father\t;remain\t;LB\t;bestLB\t;bestUB\t;status')
        self.bestUB = 57  # compact.instance.K
        self.totalNodes = 1
        self.UBs = []
        self.begin = time()
        self.timeout = 3600
        # adding optimize relaxed problem and adding cuts
        node = Node(compact.model, compact.instance)
        # compact.relax()
        # compact.model.write('teste1.lp')
        node.optNode0()
        # node.mip.write('{}_cuts.lp'.format(node.instance.instancename))
        # input('teste')
        end = time() - self.begin
        node.depth = 0
        node.node = 1
        node.father = 1
        self.eps = 0.00000001
        self.bestLB = node.mip.objective_value
        print('{:.4f}\t;{}\t;{}\t;{}\t;{:.4f};{:.4f};{:.4f}\t;{}'.format(end, node.node, '-', '0', node.mip.objective_value,
                                                              self.bestLB, self.bestUB, 'root'))
        # node.mip.write('testenode.lp')
        # input('testenode')
        self.L = collections.deque([(node.mip.objective_value, node)])


        while len(self.L) != 0:

            p = self.L.popleft()[1]  #pegar menor valor de LB da lista L
            lb = p.mip.objective_value
            if abs(lb - self.bestLB) > self.eps:
                self.bestLB = lb
            # print('Total Nodes: {} ; Remaning nodes: {}; depth: {} ; progress: LB = {} ; bestUB = {} ; gap = {}%'.
            #                                   format(self.totalNodes, len(self.L), p.depth, lb, self.bestUB,
            #                                          round((self.bestUB - lb) / self.bestUB * 100, 2)))
            set_branch_y = {}
            set_branch_set_x = {}
            set_branch_x = {}
            # input('iniciando dicionario')
            cons = [c.name for c in p.mip.constrs if c.name.startswith('depth')]
            
            for i in range(p.instance.m):
                for j in range(p.instance.n):
                    for k in range(j + 1, p.instance.n):
                        dur = time() - self.begin
                        if dur > self.timeout:
                            print(
                                'Best UB: {} ; Best LB: {}; Total nodes: {} ; Visited nodes: {} ; Remaining nodes: {} ; Integer solutions found: {}'.format(
                                    self.bestUB, self.bestLB, self.totalNodes, self.totalNodes - len(self.L), len(self.L), self.UBs
                                ))
                            print('Time: ', dur)
                            exit()
                        s = (i, j, k)
                        if isinstance(p.y[k][j][i], int) or isinstance(p.y[j][k][i], int):
                            continue
                        if (p.y[j][k][i].x + p.y[k][j][i].x) >= 0.99999999 and self.is_whole(p.y[j][k][i].x) \
                                and self.is_whole(p.y[k][j][i].x):
                            continue


                        valuejk = p.valueOrder(j, k, i)
                        valuekj = p.valueOrder(k, j, i)
                        # print('{} {} {}'.format(p.y[j][k][i].x, p.y[k][j][i].x, dist))

                        # check if is already on the model
                        # onModel = False
                        # for c in cons:
                        #     found1 = c.endswith('({},{},{})'.format(j, k, i))
                        #     found2 = c.endswith('({},{},{})'.format(k, j, i))
                        #     if found1 or found2:
                        #         onModel = True
                        #         print('exists {}'.format(c))
                        #         break
                        # if onModel:
                        #     continue

                        value = valuejk + valuekj + 5*min(valuejk, valuekj)
                        if value <= self.eps:  # only values that changes the two sides must be branched
                            # dist = min(abs(p.y[j][k][i].x), abs(p.y[k][j][i].x))
                            # # input('max(valuejk, valuekj) == 0 ; j = {} k = {} i = {} valuejk = {} valuekj = {}'.format(j, k, i, valuejk, valuekj))
                            # # input('{} {} {} {} '.format(p.y[j][k][i].x, p.y[k][j][i].x, dist, 0.5 - dist))
                            # if dist <= 0.4999:
                            #     # input()
                            #     value = dist
                            # else:
                            continue
                        # print(value)

                        # value = abs(valuejk - valuekj) / (2*max(valuejk, valuekj))  # ficar valor entre 0 e 0,5
                        # input('{} {} {} {} {} {} {}'.format(j, k, i, valuejk, valuekj, value, dist))
                        # if value > 0:
                        #     print(s)
                        # print(value, dist)
                        set_branch_y[s] = value
                        # dist = dist + value
                        # if value < 0.75:
                        #     dict[s] = dist
                        # dict[s] = value
            # print(dict)
            if len(set_branch_y) > 0:  # caso tenha variáveis y com valor diferente de 1 e 0
                set_branch_y = sorted(set_branch_y.items(), key=lambda l: l[1], reverse=True)
                # print(set_branch_y)
            else:
                for i in range(p.instance.m):
                    sum = 0
                    for j in range(p.instance.n):
                        s = (i, j)
                        # verifica se o valor de x não é inteiro. Se for, colocar no set_branch_x
                        if not self.is_whole(p.x[j][i].x):
                            set_branch_x[s] = p.x[j][i].x
                        # print('+ {} ({}) '.format(p.x[j][i].x, p.x[j][i]), end='')
                        sum += p.x[j][i].x
                    # print('= {}'.format(sum))
                    # verifica se a soma das variáveis é inteira. Caso não seja, incluir no set_branch_set_y
                    if not self.is_whole(sum):
                        s = (i, tuple(range(p.instance.n)))
                        # print('sum', sum)
                        set_branch_set_x[s] = sum
                if len(set_branch_set_x) > 0:
                    set_branch_set_x = sorted(set_branch_set_x.items(), key=lambda l: modf(l[1])[0], reverse=True)
                else:
                    set_branch_x = sorted(set_branch_x.items(), key=lambda l: abs(modf(p.x[j][i].x)[0] - 0.5), reverse=True)
            # print('ordenado')
            # print(set_branch_y, type(set_branch_y))
            # print(set_branch_set_x, type(set_branch_set_x))
            # print(set_branch_x, type(set_branch_x))
            # print('y: ', len(set_branch_y), end='')
            # print(' set_x: ', len(set_branch_set_x), end='')
            # print(' x: ', len(set_branch_x))
            # if len(set_branch_y) == 0:
            #     input()
            # input()
            self.qtd = 0
            start_iter = time()
            if len(set_branch_y) > 0:
                # print(set_branch_y)
                self.checkNodeBranchY(p, set_branch_y)
            elif len(set_branch_set_x) > 0:
                # print(set_branch_set_x)
                self.checkNodeBranchSetX(p, set_branch_set_x)
            else:
                # input()
                self.checkNodeBranchX(p, set_branch_x)
            self.L = collections.deque(sorted(self.L, key=lambda l: l[0], reverse=False))
            # print(self.L)
            end_iter = time()
            # print('Time elapsed for iteration: {}s'.format(round(end_iter - start_iter, 2)))

    def integerSol(self, p):
        for i in range(p.instance.m):
            for j in range(p.instance.n):
                if not self.is_whole(p.x[j][i].x):
                    # print(' ; var not integer: {} value {}'.format(p.x[j][i].name, p.x[j][i].x))
                    return False
                for k in range(p.instance.n):
                    if not isinstance(p.y[j][k][i], int) and not self.is_whole(p.y[j][k][i].x):
                        # print(' ; var not integer: {} value {}'.format(p.y[j][k][i].name, p.y[j][k][i].x))
                        return False
        return True

    def is_whole(self, f):
        return abs(f - round(f)) < abs(self.eps)

    def checkNodeBranchY(self, p, set_branch_y):
        self.eps = 0.00000001

        for key, value in set_branch_y:
            dur = time() - self.begin
            if dur > self.timeout:
                print(
                    'Best UB: {} ; Total nodes: {} ; Visited nodes: {} ; Remaining nodes: {} ; Integer solutions found: {}'.format(
                        self.bestUB, self.totalNodes, self.totalNodes - len(self.L), len(self.L), self.UBs
                    ))
                print('Time: ', dur)
                exit()
            i = key[0]
            j = key[1]
            k = key[2]
            # print(key, value)
            self.qtd = self.qtd + 1
            p1 = Node(p.mip, p.instance)  # j before k in machine i
            p1.depth = p.depth + 1
            p1.totalCG = p.totalCG
            p1.node = self.totalNodes + 1
            self.totalNodes += 1
            p1.father = p.node
            # input('atualizando nós')
            start = time()
            p1.updateNodeBranchY(j, k, i)
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p1.node, p1.father), end='')
            # print('Time elapsed relaxation p1: {}s ; '.format(round(end - start, 2)), end='')
            # input('teste update node')
            # p1.printSolution()

            p1.mip.write('node{}.lp'.format(p1.node))
            self.checkNode(p1)

            p2 = Node(p.mip, p.instance)  # k before j in machine i
            p2.depth = p.depth + 1
            p2.node = self.totalNodes + 1
            p2.totalCG = p.totalCG
            self.totalNodes += 1
            p2.father = p.node
            start = time()
            p2.updateNodeBranchY(k, j, i)
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p2.node,  p2.father), end='')
            # print('Time elapsed relaxation p2: {}s ; '.format(round(end - start, 2)), end='')

            p2.mip.write('node{}.lp'.format(p2.node))
            self.checkNode(p2)

            # print('quantidade de nós vistos no depth {}: {} de {}'.format(p.depth, self.qtd, len(set_branch_y)))
            # print('list size: {}'.format(len(self.L)))

    def checkNode(self, p):
        self.eps = 0.00000001
        if p.mip.status == OptimizationStatus.OPTIMAL or p.mip.status == OptimizationStatus.FEASIBLE:  # optimal or feasible solution found
            # print('checking node {}: obj p: {}'.format(p.node, p.mip.objective_value), end='')
            if p.mip.objective_value - self.bestUB <= self.eps:
                if self.is_whole(p.mip.objective_value):
                    # print(' obj integer', end='')
                    if self.integerSol(p):
                        # verificar se é melhor e então atualizar
                        # print(' ; solução inteira encontrada: {}'.format(p.mip.objective_value))
                        if p.mip.objective_value - (self.bestUB - 1) < self.eps:
                            self.bestUB = p.mip.objective_value
                            # remover de L os que estão violados (bound >= bestUB)
                            self.cutNodes()
                            p.mip.write('{}_node{}.sol'.format(p.instance.instancename.translate(str.maketrans('', '', string.punctuation)),
                                                               p.node ))
                            p.mip.write('{}_node{}.lp'.format(p.instance.instancename.translate(str.maketrans('', '', string.punctuation)),
                                                               p.node ))
                            print('{}\t;{:.4f};{:.4f};{:.4f};best integer'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
                        else:
                            print('{}\t;{:.4f};{:.4f};{:.4f};integer'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
                        self.UBs.append(p.mip.objective_value)
                        # p.printSolution()
                        # input(' check')
                    else:
                        if p.mip.objective_value - (self.bestUB - 1) < self.eps:
                            # self.totalNodes += 1
                            # print(p.mip.objective_value, self.bestUB, p.mip.objective_value - self.bestUB - 1, self.eps)
                            self.L.append((p.mip.objective_value, p))
                            print('{}\t;{:.4f};{:.4f};{:.4f};branch'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
                        else:
                            print(
                                '{}\t;{:.4f};{:.4f};{:.4f};cut'.format(len(self.L), p.mip.objective_value, self.bestLB,
                                                                       self.bestUB))
                else:
                    if p.mip.objective_value - self.bestUB < self.eps:
                        # self.totalNodes += 1
                        self.L.append((p.mip.objective_value, p))
                        print('{}\t;{:.4f};{:.4f};{:.4f};branch'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
                    else:
                        print('{}\t;{:.4f};{:.4f};{:.4f};cut'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
            else:
                print('{}\t;{:.4f};{:.4f};{:.4f};cut'.format(len(self.L), p.mip.objective_value, self.bestLB, self.bestUB))
        else:
            print('{}\t;-\t;{:.4f};{:.4f};infeasible'.format(len(self.L), self.bestLB, self.bestUB))

    def checkNodeBranchSetX(self, p, set_branch_set_x):
        self.eps = 0.00000001
        for key, value in set_branch_set_x:
            dur = time() - self.begin
            if dur > self.timeout:
                print(
                    'Best UB: {} ; Total nodes: {} ; Visited nodes: {} ; Remaining nodes: {} ; Integer solutions found: {}'.format(
                        self.bestUB, self.totalNodes, self.totalNodes - len(self.L), len(self.L), self.UBs
                    ))
                print('Time: ', dur)
                exit()
            self.qtd = self.qtd + 1
            # input('atualizando nós')
            p1, p2 = p.updateNodeBranchSetX(key[0], key[1], value)
            p1.node = self.totalNodes + 1
            p1.totalCG = p.totalCG
            self.totalNodes += 1
            p1.father = p.node
            p2.node = self.totalNodes + 1
            self.totalNodes += 1
            p2.father = p.node
            p2.totalCG = p.totalCG
            start = time()
            p1.optNode()
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p1.node, p1.father), end='')

            # print('Time elapsed relaxation p1: {}s ; '.format(round(end - start, 2)), end='')
            # input('teste update node')
            # p1.printSolution()

            # p1.mip.write('node{}.lp'.format(p1.node))
            self.checkNode(p1)

            start = time()
            p2.optNode()
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p2.node, p2.father), end='')
            # print('Time elapsed relaxation p2: {}s ; '.format(round(end - start, 2)), end='')

            # p2.mip.write('node{}.lp'.format(p2.node))
            self.checkNode(p2)

            # remover de L os que estão violados (bound >= bestUB)
            # self.cutNodes()

            # print('quantidade de nós vistos no depth {}: {} de {}'.format(p.depth, self.qtd, len(set_branch_set_x)))
            # print('list size: {}'.format(len(self.L)))

    def cutNodes(self):
        self.eps = 0.00000001
        for index in range(len(self.L) - 1, -1, -1):
            if self.L[index][0] - (self.bestUB - 1) > self.eps:
                # print(index, self.L[index].mip.objective_value)
                # input()
                del self.L[index]

    def checkNodeBranchX(self, p, set_branch_x):
        self.eps = 0.00000001
        for key, value in set_branch_x:
            dur = time() - self.begin
            if dur > self.timeout:
                print(
                    'Best UB: {} ; Total nodes: {} ; Visited nodes: {} ; Remaining nodes: {} ; Integer solutions found: {}'.format(
                        self.bestUB, self.totalNodes, self.totalNodes - len(self.L), len(self.L), self.UBs
                    ))
                print('Time: ', dur)
                exit()
            self.qtd = self.qtd + 1
            # input('atualizando nós')
            # print(key, value)
            p1, p2 = p.updateNodeBranchX(key[0], key[1], value)
            p1.node = self.totalNodes + 1
            self.totalNodes += 1
            p2.node = self.totalNodes + 1
            self.totalNodes += 1
            start = time()
            p1.optNode()
            end = time()
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p1.node, p1.father), end='')
            # print('Time elapsed relaxation p1: {}s ; '.format(round(end - start, 2)), end='')
            # input('teste update node')
            # p1.printSolution()

            # p1.mip.write('testep1depth{}.lp'.format(p1.depth))
            self.checkNode(p1)

            start = time()
            p2.optNode()
            end = time()
            dur = end - self.begin
            # print('time\t;node\t;remain\t;father\t;LB\t;bestUB\t;status')
            print('{:.4f}\t;{}\t;{}\t;'.format(dur, p2.node, p2.father), end='')
            # print('Time elapsed relaxation p2: {}s ; '.format(round(end - start, 2)), end='')

            # p2.mip.write('testep2depth{}.lp'.format(p2.depth))
            self.checkNode(p2)

            # remover de L os que estão violados (bound >= bestUB)
            # self.cutNodes()

            # print('quantidade de nós vistos no depth {}: {} de {}'.format(p.depth, self.qtd, len(set_branch_x)))
            # print('list size: {}'.format(len(self.L)))
