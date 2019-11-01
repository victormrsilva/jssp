from time import time
from itertools import combinations, permutations

import numpy as np

from compact import Compact
from mip import Model, xsum, BINARY, INTEGER, CONTINUOUS, OptimizationStatus
import collections


class Node:
    def __init__(self, problem, instance):
        self.depth = 0
        # same instance always
        self.instance = instance

        # auxiliar model
        self.m_aux = Model(solver_name="cbc")
        self.m_aux.verbose = 0

        # variables for cuts
        self.iterationsCuts = 0
        self.lc = 2  # lower number of jobs in a clique to be found
        self.hc = 10  # high number of jobs in a clique to be found
        self.maxCliques = 100  # maximum number of cliques to be investigated
        if self.instance.n < 10:
            self.hc = self.instance.n
        self.totalBasicCuts = 0
        self.totalBasicCutsEpsilon = 0
        self.totalTriangleCuts = 0
        self.totalCliqueCuts = 0
        self.totalTwoJobsCuts = 0
        self.totalLateJobCuts = 0
        self.totalHalfCuts = 0

        # copy the model
        self.mip = problem.copy()

        # model variables
        self.c = 0
        self.x = [[0 for i in range(self.instance.m)] for j in range(self.instance.n)]
        self.y = [
            [[0 for i in range(self.instance.m)] for k in range(self.instance.n)] for j in range(self.instance.n)]
        self.v = [
            [[0 for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]
        self.u = [
            [[0 for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]

        for i, var in enumerate(self.mip.vars):
            if var.name.startswith('x'):
                pos = var.name.split('(')[1].split(')')[0].split(',')  # get j, i
                self.x[int(pos[0])][int(pos[1])] = var
            elif var.name.startswith('y'):
                pos = var.name.split('(')[1].split(')')[0].split(',')  # get j, k, i
                self.y[int(pos[0])][int(pos[1])][int(pos[2])] = var
            elif var.name.startswith('C'):
                self.c = var
            elif var.name.startswith('v'):
                pos = var.name.split('(')[1].split(')')[0].split(',')  # get j, k, i
                self.v[int(pos[0])][int(pos[1])][int(pos[2])] = var
            elif var.name.startswith('u'):
                pos = var.name.split('(')[1].split(')')[0].split(',')  # get j, k, i
                self.u[int(pos[0])][int(pos[1])][int(pos[2])] = var

    # run optimizarions until all cuts have been found
    def optNode0(self):
        self.mip.relax()
        start = time()

        self.mip.verbose = 0
        self.mip.optimize()
        # self.printSolution()

        # check for clique cuts
        for a in range(self.instance.m):
            # clique_cuts = 0
            clique_cuts = self.clique_cuts_best(a)
            self.totalCliqueCuts += clique_cuts
        newConstraints = True

        while newConstraints:
            self.iterationsCuts += 1
            newConstraints = False
            lastobj = self.mip.objective_value
            # start = time()
            hasCuts = self.cuts()
            if hasCuts > 0:
                newConstraints = True
                self.mip.optimize()

        end = time()
        # print('Time elapsed relaxation: {}s'.format(round(end - start, 2)))


    # j before k in machine i
    def updateNode(self, j, k, i):
        boundchange = False
        if self.x[j][i].lb >= self.x[k][i].lb:
            boundchange = True
            # print('est atualizado')
            self.x[k][i].lb = self.x[j][i].lb + self.instance.times[j][i]
            # atualizar das maquinas seguintes
            for a in range(self.instance.o[k][i] + 1, self.instance.m):
                o = self.instance.machines[k][a]
                o_a = self.instance.machines[k][a - 1]
                self.x[k][o].lb = self.x[k][o_a].lb + self.instance.times[k][o_a]

        if self.x[j][i].ub >= self.x[k][i].ub:
            boundchange = True
            self.x[j][i].ub = self.x[k][i].ub - self.instance.times[j][i]
            ub = self.x[k][i].ub - self.instance.times[j][i]
            # atualizar das maquinas anteriores
            for a in range(self.instance.o[j][i] - 1, -1, -1):
                o = self.instance.machines[j][a]
                ub = ub - self.instance.times[j][o]
                self.x[j][o].ub = ub

        self.mip += self.x[k][i] - self.x[j][i] >= self.instance.times[j][i], 'depth{}({},{},{})'.format(self.depth, j, k, i)

        # caso tenha mudança do bound, atualizar mccormick
        if boundchange:
            self.updateMcCormick(j, k, i)
        # atualiza nó com novos cortes
        self.optNode()

    def valueOrder(self, j, k, i):
        value1 = 0
        if self.x[j][i].lb >= self.x[k][i].lb:
            # print('est atualizado')
            lb_oa = self.x[j][i].lb + self.instance.times[j][i]
            value1 += lb_oa - self.x[k][i].lb
            # atualizar das maquinas seguintes
            for a in range(self.instance.o[k][i] + 1, self.instance.m):
                o = self.instance.machines[k][a]
                o_a = self.instance.machines[k][a - 1]
                lb_oa = lb_oa + self.instance.times[k][o_a]
                value1 += lb_oa - self.x[k][o].lb
        value2 = 0
        if self.x[j][i].ub >= self.x[k][i].ub:
            ub_od = self.x[k][i].ub - self.instance.times[j][i]  # self.x[j][i].ub
            # atualizar das maquinas anteriores
            for a in range(self.instance.o[j][i] - 1, -1, -1):
                o = self.instance.machines[j][a]
                ub_od = ub_od - self.instance.times[j][o]
                # print('o = {}, o_d = {}'.format(o, o_d))
                value2 += self.x[j][o].ub - ub_od

        value = value1 + value2
        # value = float(value / (2 * self.instance.K))
        return value

    def updateMcCormick(self, j, k, i):
        if k < j:
            k, j = j, k

        for i in range(self.instance.m):
            vjkL = self.x[k][i].lb - self.x[j][i].ub - self.instance.times[j][i]
            vjkU = self.x[k][i].ub - self.x[j][i].lb - self.instance.times[j][i]
            vkjL = self.x[j][i].lb - self.x[k][i].ub - self.instance.times[k][i]
            vkjU = self.x[j][i].ub - self.x[k][i].lb - self.instance.times[k][i]
            names_constrs = ['MCLin1({},{},{})'.format(j, k, i),
                             'MCLin2({},{},{})'.format(j, k, i),
                             'MCLin3({},{},{})'.format(j, k, i),
                             'MCLin4({},{},{})'.format(j, k, i),
                             'MCLin5({},{},{})'.format(j, k, i),
                             'MCLin6({},{},{})'.format(k, j, i),
                             'MCLin7({},{},{})'.format(k, j, i),
                             'MCLin8({},{},{})'.format(k, j, i),
                             'MCLin9({},{},{})'.format(k, j, i),
                             'MCLin10({},{},{})'.format(k, j, i)]
            # remover constraints antigas
            self.mip.remove([c for c in self.mip.constrs if c.name in names_constrs])

            # refazer constraints removidas
            self.mip += self.u[j][k][i] >= 0, 'MCLin1({},{},{})'.format(j, k, i)
            self.mip += self.u[j][k][i] - vjkL * self.y[j][k][i] >= 0, 'MCLin2({},{},{})'.format(j, k, i)
            self.mip += self.u[j][k][i] - vjkU * self.y[j][k][i] - self.v[j][k][i] >= - vjkU, 'MCLin3({},{},{})'.format(j, k, i)
            self.mip += self.u[j][k][i] - vjkU * self.y[j][k][i] <= 0, 'MCLin4({},{},{})'.format(j, k, i)
            self.mip += self.u[j][k][i] - self.v[j][k][i] - vjkL * self.y[j][k][i] <= - vjkL, 'MCLin5({},{},{})'.format(j, k, i)
            self.mip += self.u[k][j][i] >= 0, 'MCLin6({},{},{})'.format(k, j, i)
            self.mip += self.u[k][j][i] + vkjL * self.y[k][j][i] >= vkjL, 'MCLin7({},{},{})'.format(k, j, i)
            self.mip += self.u[k][j][i] + vkjU * self.y[j][k][i] - self.v[k][j][i] >= 0, 'MCLin8({},{},{})'.format(k, j, i)
            self.mip += self.u[k][j][i] - vkjU * self.y[j][k][i] <= vkjU, 'MCLin9({},{},{})'.format(k, j, i)
            self.mip += self.u[k][j][i] + vkjL * self.y[j][k][i] - self.v[k][j][i] <= 0, 'MCLin10({},{},{})'.format(k, j, i)

    def optNode(self):
        self.mip.relax()
        start = time()

        self.mip.verbose = 0
        self.mip.optimize()
        # self.printSolution()
        # if infesiable return
        # print(self.mip.status)
        # self.mip.write('teste.lp')
        if self.mip.status != OptimizationStatus.OPTIMAL and self.mip.status != OptimizationStatus.FEASIBLE:
            return
        # remove constraints with slack >= 0.5
        setconstr = [c for c in self.mip.constrs if c.name.startswith('cut') and c.slack >= abs(0.5)]
        self.mip.remove(setconstr)
        # for c in setconstr:
        #     self.checkConstraints(c)
        #     input()
        #     print('slack: {} ; constr: {}'.format(c.slack, c))
        # input('constraints')

        # check for clique cuts
        for a in range(self.instance.m):
            # clique_cuts = 0
            clique_cuts = self.clique_cuts_best(a)
            self.totalCliqueCuts += clique_cuts

        newConstraints = True
        while newConstraints:
            self.iterationsCuts += 1
            newConstraints = False
            lastobj = self.mip.objective_value
            # start = time()
            hasCuts = self.cuts()
            if hasCuts > 0:
                newConstraints = True
                self.mip.optimize()
            # self.cuts()
            #
            # self.mip.optimize()
            # # print('lastobj = {} ; obj = {}'.format(lastobj, self.mip.objective_value))
            #
            # # apenas verifica novamente se há mudança no objetivo
            # if (self.mip.objective_value - lastobj) >= 0.000001:
            #     newConstraints = True

        end = time()
        # print('Time elapsed relaxation: {}s'.format(round(end - start, 2)))

    def cuts(self):
        totalCuts = 0
        for a in range(self.instance.m):
            # triangle_cuts = 0
            triangle_cuts = self.triangle_cuts_best(a)

            # basic_cuts = 0
            basic_cuts = self.basic_cuts_best(a)

            # two_job_cuts = 0
            two_job_cuts = self.two_jobs_cuts_best(a)

            late_jobs_cuts = 0
            half_cuts = 0
            basic_cuts_epsilon = 0
            for k in range(self.instance.n):
                half_cuts += self.half_cuts_best(a, k)

                basic_cuts_epsilon = self.basic_cuts_plus_epsilon_best(a, k)

                for l in range(self.instance.n):
                    late_jobs_cuts += self.late_job_cuts_best(a, k, l)

            self.totalLateJobCuts += triangle_cuts
            self.totalBasicCuts += basic_cuts
            self.totalTwoJobsCuts += two_job_cuts
            self.totalHalfCuts += half_cuts
            self.totalBasicCutsEpsilon += basic_cuts_epsilon
            self.totalLateJobCuts += late_jobs_cuts
            totalCuts += basic_cuts + basic_cuts_epsilon + half_cuts + late_jobs_cuts + two_job_cuts + triangle_cuts
        return totalCuts

    def checkConstraints(self, c):
        print(c)
        soma = 0
        line = ''
        for (var, val) in c.expr.expr.items():
            astr = ' {:+}*{}'.format(val, var.x)
            line += astr
            soma += val * var.x

        rhs = c.expr.const * -1.0
        dir = ''
        if c.expr.sense == '=':
            dir += ' = {}'.format(rhs)
        elif c.expr.sense == '<':
            dir += ' <= {}'.format(rhs)
        elif c.expr.sense == '>':
            dir += ' >= {}'.format(rhs)
        print('{} {}'.format(line, dir))

        print('{} {}. Slack: {}'.format(soma, dir, c.slack))

    def E(self, S, a):
        min = 99999999
        for s in S:
            if self.instance.e[s][a] < min:
                min = self.instance.e[s][a]
        return min

    def F(self, S, a):
        min = 99999999
        for s in S:
            if self.instance.f[s][a] < min:
                min = self.instance.f[s][a]
        return min

    def p(self, S, a):
        soma = 0
        for s in S:
            soma += self.instance.times[s][a]
        return soma

    def basic_cuts_best(self, a):
        m = self.m_aux
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        e_p = [m.add_var(var_type=INTEGER, lb=0, name='e_p({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        y = [m.add_var(var_type=INTEGER, lb=0, name='y({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='e({})'.format(i)) for i in range(self.instance.n)]
        z = m.add_var(var_type=INTEGER, name='E')

        m.objective = xsum(
            self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n)) - xsum(
            self.instance.times[j][a] * e_p[j] for j in range(self.instance.n)) - xsum(
            self.instance.times[j][a] * self.instance.times[i][a] * xij[j][i] for i in range(self.instance.n) for j in
            range(i + 1, self.instance.n))

        for j in range(self.instance.n):
            m += v[j] - self.instance.e[j][a] * x_aux[j] + self.instance.K * x_aux[
                j] == self.instance.K, 'eq26({})'.format(j)
            m += z - v[j] <= 0, 'eq27({})'.format(j)
            m += z - y[j] >= 0, 'eq28({})'.format(j)
            m += y[j] - self.instance.K * o[j] <= 0, 'eq29({})'.format(j)
            m += y[j] - v[j] <= 0, 'eq30({})'.format(j)
            m += y[j] - v[j] - self.instance.K * o[j] >= - self.instance.K, 'eq31({})'.format(j)
            # z*x(j)
            m += e_p[j] - self.instance.K * x_aux[j] <= 0, 'e_s1({})'.format(j)
            m += e_p[j] - z <= 0, 'e_s2({})'.format(j)
            m += e_p[j] - z - self.instance.K * x_aux[j] >= -self.instance.K, 'e_s3({})'.format(j)
            # x[i]*x[j]
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += xsum(o[j] for j in range(self.instance.n)) == 1, 'eq32'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        # left side
        ls = self.E(S, a) * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        c_name = 'cut_basic_best({},{})'.format(''.join(str(i) for i in S), a)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)
        self.mip += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, c_name

        # self.mip += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, 'cut_basic_best' \
        #                                                                              '{}({},{})'.format(
        #     self.iterationsCuts, ''.join(str(i) for i in S), a)
        return 1

    def basic_cuts_plus_epsilon_best(self, a, k):
        m = self.m_aux
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                xij[i][j] = xij[j][i]

        var = xsum(
            self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n))  # part 1 of cut
        var += - xsum(
            self.instance.times[j][a] * self.instance.times[i][a] * xij[i][j] for i in range(self.instance.n) for j in
            range(i + 1, self.instance.n))  # part 2 of cut
        var += - self.instance.e[k][a] * xsum(
            self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k)  # part 3 of cut
        var += xsum(
            self.y[j][k][a].x * max(self.instance.e[k][a] - self.instance.e[j][a], 0) * self.instance.times[j][a] *
            x_aux[j] for j in range(self.instance.n) if
            j != k)  # part 4 of cut if i == j (just expand the multiplication)
        var += xsum(
            self.y[i][k][a].x * max(self.instance.e[k][a] - self.instance.e[i][a], 0) * self.instance.times[j][a] *
            xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if
            i != k and j != k and i != j)  # part 4 of cut if i != j (just expand the multiplication)
        m.objective = var

        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += x_aux[k] == 0, 'not_in_S'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        # left side
        ls = self.instance.e[k][a] * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        c_name = 'cut_basic_plus_epsilon_best({},{},{})'.format(''.join(str(i) for i in S), a, k)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)
        self.mip += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
            self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,
                                                                                                         a) >= ls, c_name

        # self.mip += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
        #     self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,
        #                                                                                                  a) >= ls, 'cut_basic_plus_epsilon_best{}({},{},{})'.format(
        #     self.iterationsCuts, ''.join(str(i) for i in S), a, k)
        return 1

    def half_cuts_best(self, a, k):
        m = self.m_aux
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        t = [m.add_var(var_type=INTEGER, lb=0, name='t({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='o({})'.format(i)) for i in range(self.instance.n)]
        e = m.add_var(var_type=INTEGER, name='E')
        C = m.add_var(name='C', lb=-100 * self.instance.K)

        m.objective = C

        for j in range(self.instance.n):
            m += v[j] - self.instance.e[j][a] * x_aux[j] + self.instance.K * x_aux[
                j] == self.instance.K, 'eq26({})'.format(j)
            m += e - v[j] <= 0, 'eq27({})'.format(j)
            m += e - t[j] >= 0, 'eq28({})'.format(j)
            m += t[j] - self.instance.K * o[j] <= 0, 'eq29({})'.format(j)
            m += t[j] - v[j] <= 0, 'eq30({})'.format(j)
            m += t[j] - v[j] - self.instance.K * o[j] >= - self.instance.K, 'eq31({})'.format(j)
        m += xsum(o[j] for j in range(self.instance.n)) == 1, 'eq32'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'
        m += C + e + xsum(
            self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k) == \
             self.x[k][a].x

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        c_name = 'cut_half_best({},{},{})'.format(k, ''.join(str(i) for i in S), a)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)
        self.mip += self.x[k][a] - self.E(S, a) - xsum(
            self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= 0, c_name

        # self.mip += self.x[k][a] - self.E(S, a) - xsum(
        #     self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= 0, 'cut_half_best{}({},{},{})'.format(
        #     self.iterationsCuts, k, ''.join(str(i) for i in S), a)
        # input()
        return 1

    def late_job_cuts_best(self, a, k, l):
        m = self.m_aux
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        C = m.add_var(name='C', lb=-100 * self.instance.K)

        m.objective = C

        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'
        m += C + xsum(self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n)
                      if j != k) - xsum(self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0)
                                        * x_aux[j] for j in range(self.instance.n) if j != l) == self.x[k][a].x - \
             self.instance.e[l][a]

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        c_name = 'cut_late_job_best({},{},{},{})'.format(''.join(str(i) for i in S), a, k, l)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)
        self.mip += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
            self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
                    self.instance.e[l][a], c_name

        # self.mip += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
        #     self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
        #             self.instance.e[l][a], 'cut_late_job_best{}({},{},{},{})'.format(self.iterationsCuts,
        #                                                                              ''.join(str(i) for i in S), a,
        #                                                                              k, l)
        return 1

    def two_jobs_cuts_best(self, a):
        # print('Best two jobs cuts')
        m = self.m_aux
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                xij[i][j] = xij[j][i]

        var = xsum(self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n))
        var += xsum(self.instance.e[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n))
        for i in range(self.instance.n):
            for j in range(self.instance.n):
                if i != j:
                    var -= self.instance.e[i][a] * self.x[j][a].x * xij[i][j]
                    var -= self.instance.e[i][a] * self.instance.times[j][a] * xij[i][j]
                if j > i:
                    var -= self.instance.times[i][a] * self.instance.times[j][a] * xij[i][j]

        m.objective = var

        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += xsum(x_aux[j] for j in range(self.instance.n)) == 2, 'must_be_2'

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        self.two_jobs_cuts(S[0], S[1], a)
        # input()
        return 1

    def two_jobs_cuts(self, i, j, a):
        # right side
        rs = (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a].x + (
                self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a].x
        # left side
        ls = self.instance.times[i][a] * self.instance.times[j][a] + self.instance.e[i][a] * \
             self.instance.times[j][a] + self.instance.e[j][a] * self.instance.times[i][a]

        # apply cut to the model
        c_name = 'cut_two_jobs({},{},{})'.format(i, j, a)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)
        self.mip += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
            a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
                        a] >= ls, c_name

        # self.mip += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
        #     a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
        #                 a] >= ls, 'cut_two_jobs{}({},{},{})'.format(self.iterationsCuts, i, j, a)
        return 1

    def clique_cuts(self, S, a):
        s_aux = list(range(len(S)))
        perms = np.asarray(list(permutations(s_aux)))
        perms = [np.asarray(s) for s in perms]

        K = [[0] * len(S) for i in range(len(perms))]

        for i in range(len(K)):
            soma = 0
            for j in range(len(K[i])):
                if j == 0:
                    soma += self.instance.e[S[perms[i][j]]][a]
                else:
                    soma += self.instance.times[S[perms[i][j - 1]]][a]
                for aux in range(j, len(K[i])):
                    K[i][perms[i][aux]] = soma
        m = self.m_aux
        m.clear()
        m.verbose = 0
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        m.objective = xsum(i for i in t)

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)

        m.optimize()

        soma = 0
        for i in range(len(S)):
            soma += t[i].x * self.x[S[i]][a].x

        cuts = 0
        # xt >= is a valid inequality. Thus, a violated cut is < 1
        if (1 - soma) > 0.00001:  # < 1:
            cuts = 1
            c_name = 'cut_clique({},{})'.format(''.join(str(i) for i in S), a)
            c = self.mip.constr_by_name(c_name)
            if c is not None:
                self.mip.remove(c)
            self.mip += xsum(t[i].x * self.x[S[i]][a] for i in range(len(S))) >= 1, c_name

        return cuts

    """
    # sc = size of smaller clique
    # lc = size of larger clique
    # maxCliques = maximum cliques to be investigated
    """

    def clique_cuts_best(self, a):
        cliquesFound = 0
        start = time()
        timeLimit = 600  # 10 minutes
        for sizeS in range(self.lc, self.hc):
            comb = combinations(list(range(0, self.instance.n)),
                                sizeS + 1)  # combinations of all possibles jobs of size sizeS
            comb = list(comb)
            dict = {}
            for s in range(len(comb)):
                end = time()
                if (end - start) >= timeLimit:
                    return 0
                dist = 0
                S = comb[s]
                for i in list(range(0, sizeS + 1)):
                    for j in list(range(i + 1, sizeS + 1)):
                        dist += abs(self.x[S[i]][a].x - self.x[S[j]][a].x)
                dict[s] = dist
            i = 0
            end = time()
            if (end - start) >= timeLimit:
                print('Time limit for enumaration of cliques')
                return 0

            for key, value in sorted(dict.items(), key=lambda l: l[1]):
                end = time()
                if (end - start) >= timeLimit:
                    print('Time limit for finding cliques. Cliques found: {}'.format(cliquesFound))
                    return cliquesFound
                S = comb[key]
                cliquesFound += self.clique_cuts(S, a)
                i += 1
                if i == self.maxCliques:
                    break
        return cliquesFound

    def triangle_cuts_best(self, a):
        # print('Best triangle cuts')
        m = self.m_aux
        m.clear()
        m.verbose = 0
        xij = [[m.add_var(var_type=BINARY, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)]
               for j in range(self.instance.n)]

        var = - xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if i != j)
        m.objective = var

        for i in range(self.instance.n):
            m += xsum(xij[i][j] for j in range(self.instance.n) if j != i) - xsum(xij[j][i]
                                                                                  for j in range(self.instance.n) if
                                                                                  j != i) == 0, 'flow({})'.format(i)
        m += xsum(
            self.y[i][j][a].x * xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i
        ) - xsum(
            xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i) >= -1, 'triangle_cut'
        m.optimize()

        if m.objective_value > -0.000000001 or m.status == OptimizationStatus.INFEASIBLE:
            return 0

        S = []
        names = []
        soma = 0
        for i in range(self.instance.n):
            for j in range(self.instance.n):
                if xij[i][j].x > 0.99999999:
                    names.append('({}{})'.format(i,j))
                    S.append(self.y[i][j][a])
                    soma += self.y[i][j][a].x * xij[i][j].x

        if len(S) <= 1 or soma - (len(S) - 1) <= 0.000000001:
            return 0
        c_name = 'cut_triangle_best({},{})'.format(''.join(str(i) for i in names), a)
        c = self.mip.constr_by_name(c_name)
        if c is not None:
            self.mip.remove(c)

        self.mip += xsum(i for i in S) <= len(S) - 1, c_name
        return 1

    def printSolution(self):
        print("C: ", self.c.x)
        for j in range(self.instance.n):
            for i in range(self.instance.m):
                print('x({},{}) = {} '.format(j, i, self.x[j][i].x), end='')
            print()
