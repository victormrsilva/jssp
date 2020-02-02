# compact.py
import string

from mip.model import Model, xsum, Var
from mip.constants import INTEGER, BINARY, CONTINUOUS, OptimizationStatus
import sys
from time import process_time, time
import numpy as np
import copy

# from compact_cutpool import Compact_CutPool
from SA_clique import Clique

from itertools import permutations, combinations


class Compact:
    def __init__(self, instance):
        self.instance = instance
        self.m = Model(solver_name="cbc")
        self.m.verbose = 0
        self.model = Model(name='compact', solver_name="cbc")
        self.model.verbose = 0
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


    # build the problem with big-M
    def constructProblemM(self):
        # self.instance.print()

        self.c = self.model.add_var(var_type=INTEGER, name="C")
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                self.x[j][i] = self.model.add_var(var_type=INTEGER, name='x(j{},m{})'.format(j, i))
                for k in range(j+1, self.instance.n):
                    self.y[j][k][i] = self.model.add_var(var_type=BINARY, name='y(j{},k{},m{})'.format(j, k, i))
                    self.y[k][j][i] = self.model.add_var(var_type=BINARY, name='y(j{},k{},m{})'.format(k, j, i))

        self.model.objective = self.c

        # constraints (2)
        for j in range(self.instance.n):
            for i in range(1, self.instance.m):
                self.model += self.x[j][self.instance.machines[j][i]] - self.x[j][self.instance.machines[j][i - 1]] >= \
                              self.instance.times[j][self.instance.machines[j][i - 1]], 'ord({},{})'.format(j, i)

        # constraints (3-4)
        for j in range(self.instance.n):
            for k in range(self.instance.n):
                if k != j:
                    for i in range(self.instance.m):
                        self.model += self.x[j][i] - self.x[k][i] + self.instance.K * self.y[j][k][i] >= \
                                      self.instance.times[k][i], 'phi({},{},{})'.format(j, k, i)
                        self.model += -self.x[j][i] + self.x[k][i] - self.instance.K * self.y[j][k][i] >= \
                                      self.instance.times[j][i] - self.instance.K, 'psy({},{},{})'.format(j, k, i)

        # constraints (5)
        for j in range(self.instance.n):
            self.model += self.c - self.x[j][self.instance.machines[j][self.instance.m - 1]] >= self.instance.times[j][
                self.instance.machines[j][self.instance.m - 1]], 'makespan({})'.format(j)
        self.model.write(
            '{}_modelM.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))

    # build the problem with McCormick linearization
    def constructProblemMcCormick(self):
        print('McCormic linearization')
        self.model.clear()
        self.c = self.model.add_var(var_type=INTEGER, lb=0, ub=self.instance.K, name="C")

        for i in range(self.instance.m):
            for j in range(self.instance.n):
                self.x[j][i] = self.model.add_var(var_type=INTEGER, lb=self.instance.est[j][i],
                                                  ub=self.instance.lst[j][i],
                                                  name='x({},{})'.format(j, i))
                for k in range(j+1, self.instance.n):
                    self.y[j][k][i] = self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(j, k, i))
                    self.y[k][j][i] = self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(k, j, i))
                    self.v[j][k][i] = self.model.add_var(var_type=INTEGER,
                                        lb=self.instance.est[j][i] - self.instance.lst[k][i] - self.instance.times[k][i],
                                        ub=self.instance.lst[j][i] - self.instance.est[k][i] - self.instance.times[k][i],
                                        name='v({},{},{})'.format(j, k, i))
                    self.v[k][j][i] = self.model.add_var(var_type=INTEGER,
                                        lb=self.instance.est[k][i] - self.instance.lst[j][i] - self.instance.times[j][i],
                                        ub=self.instance.lst[k][i] - self.instance.est[j][i] - self.instance.times[j][i],
                                        name='v({},{},{})'.format(k, j, i))
                    self.u[j][k][i] = self.model.add_var(var_type=INTEGER,
                                                         lb=-10*self.instance.K,
                                                         name='u({},{},{})'.format(j, k, i))
                    self.u[k][j][i] = self.model.add_var(var_type=INTEGER,
                                                         lb=-10*self.instance.K,
                                                         name='u({},{},{})'.format(k, j, i))


        self.model.objective = self.c

        # constraints (2)
        for j in range(self.instance.n):
            for i in range(1, self.instance.m):
                self.model += self.x[j][self.instance.machines[j][i]] - self.x[j][self.instance.machines[j][i - 1]] >= \
                              self.instance.times[j][self.instance.machines[j][i - 1]], 'ord({},{})'.format(j, i)

        for j in range(self.instance.n):
            for k in range(j + 1, self.instance.n):
                for i in range(self.instance.m):
                    self.model += self.y[j][k][i] + self.y[k][j][i] == 1, 'triangle2({},{},{})'.format(j, k, i)


        # constraints (3-4)
        for j in range(self.instance.n):
            for k in range(j + 1, self.instance.n):
                for i in range(self.instance.m):
                    vjkL = self.v[j][k][i].lb
                    vjkU = self.v[j][k][i].ub
                    vkjL = self.v[k][j][i].lb
                    vkjU = self.v[k][j][i].ub
                    self.model += self.v[j][k][i] == self.x[j][i] - self.x[k][i] - self.instance.times[k][i], 'MCLinV({},{},{})'.format(j, k, i)
                    self.model += self.v[k][j][i] == self.x[k][i] - self.x[j][i] - self.instance.times[j][i], 'MCLinV({},{},{})'.format(k, j, i)
                    self.model += self.u[j][k][i] >= 0, 'MCLin1({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] >= vjkL * self.y[k][j][i], 'MCLin2({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] >= self.v[j][k][i] + vjkU * self.y[k][j][i] - vjkU, 'MCLin3({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] <= vjkU * self.y[k][j][i], 'MCLin4({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] <= self.v[j][k][i] + vjkL * self.y[k][j][i] - vjkL, 'MCLin5({},{},{})'.format(j, k, i)
                    self.model += self.u[k][j][i] >= 0, 'MCLin6({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] >= vkjL * self.y[j][k][i], 'MCLin7({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] >= self.v[k][j][i] + vkjU * self.y[j][k][i] - vkjU, 'MCLin8({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] <= vkjU * self.y[j][k][i], 'MCLin9({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] <= self.v[k][j][i] + vkjL * self.y[j][k][i] - vkjL, 'MCLin10({},{},{})'.format(k, j, i)

        # constraints (5)
        for j in range(self.instance.n):
            self.model += self.c - self.x[j][self.instance.machines[j][self.instance.m - 1]] >= self.instance.times[j][
                self.instance.machines[j][self.instance.m - 1]], 'makespan({})'.format(j)
        self.model.write(
            '{}_modelMCLin.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))

    def constructProblemMcCormickNonNegative(self):
        # print('McCormic linearization non-negative')

        self.c = self.model.add_var(var_type=INTEGER, lb=0, ub=self.instance.K, name="C")

        for i in range(self.instance.m):
            for j in range(self.instance.n):
                self.x[j][i] = self.model.add_var(var_type=INTEGER, lb=self.instance.est[j][i],
                                                  ub=self.instance.lst[j][i],
                                                  name='x({},{})'.format(j, i))
                for k in range(j+1, self.instance.n):
                    self.y[j][k][i] = self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(j, k, i))
                    self.y[k][j][i] = self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(k, j, i))

                    lb = self.instance.est[j][i] - self.instance.lst[k][i] - self.instance.times[k][i]
                    ub = self.instance.lst[j][i] - self.instance.est[k][i] - self.instance.times[k][i]
                    if lb >= 0:
                        self.v[j][k][i] = (self.model.add_var(var_type=INTEGER,
                                        lb=lb,
                                        ub=ub,
                                        name='v({},{},{})p'.format(j, k, i)), 0)
                    else:
                        self.v[j][k][i] = (self.model.add_var(var_type=INTEGER,
                                        lb=0,
                                        ub=ub,
                                        name='v({},{},{})p'.format(j, k, i)),
                                            self.model.add_var(var_type=INTEGER,
                                        lb=0,
                                        ub=abs(lb),
                                        name='v({},{},{})n'.format(j, k, i)))
                    lb = self.instance.est[k][i] - self.instance.lst[j][i] - self.instance.times[j][i]
                    ub = self.instance.lst[k][i] - self.instance.est[j][i] - self.instance.times[j][i]
                    if lb >= 0:
                        self.v[k][j][i] = (self.model.add_var(var_type=INTEGER,
                                        lb=lb,
                                        ub=ub,
                                        name='v({},{},{})p'.format(k, j, i)), 0)
                    else:
                        self.v[k][j][i] = (self.model.add_var(var_type=INTEGER,
                                            lb=0,
                                            ub=ub,
                                            name='v({},{},{})p'.format(k, j, i)),
                                           self.model.add_var(var_type=INTEGER,
                                            lb=0,
                                            ub=abs(lb),
                                            name='v({},{},{})n'.format(k, j, i))
                                          )
                    self.u[j][k][i] = (self.model.add_var(var_type=INTEGER,
                                                         lb=0,
                                                         ub=10*self.instance.K,
                                                         name='u({},{},{})p'.format(j, k, i)),
                                       self.model.add_var(var_type=INTEGER,
                                                          lb=0,
                                                          ub=10 * self.instance.K,
                                                          name='u({},{},{})n'.format(j, k, i))
                                       )
                    self.u[k][j][i] = (self.model.add_var(var_type=INTEGER,
                                                         lb=0,
                                                         ub=10*self.instance.K,
                                                         name='u({},{},{})p'.format(k, j, i)),
                                       self.model.add_var(var_type=INTEGER,
                                                          lb=0,
                                                          ub=10 * self.instance.K,
                                                          name='u({},{},{})n'.format(k, j, i))
                                       )

        self.model.objective = self.c

        # constraints (2)
        for j in range(self.instance.n):
            for i in range(1, self.instance.m):
                self.model += self.x[j][self.instance.machines[j][i]] - self.x[j][self.instance.machines[j][i - 1]] >= \
                              self.instance.times[j][self.instance.machines[j][i - 1]], 'ord({},{})'.format(j, i)

        for j in range(self.instance.n):
            for k in range(j + 1, self.instance.n):
                for i in range(self.instance.m):
                    self.model += self.y[j][k][i] + self.y[k][j][i] == 1, 'triangle2({},{},{})'.format(j, k, i)


        # constraints (3-4)
        for j in range(self.instance.n):
            for k in range(j + 1, self.instance.n):
                for i in range(self.instance.m):
                    if type(self.v[j][k][i][1]) == int:  # indica que a variável não possui lb negativo
                        vjkL = self.v[j][k][i].lb
                        vjkU = self.v[j][k][i].ub
                    else: # indica que a variável possui lb negativo
                        vjkL = -self.v[j][k][i][1].ub  # lb
                        vjkU = self.v[j][k][i][0].ub  # ub

                    if type(self.v[k][j][i][1]) == int: # indica que a variável não possui lb negativo
                        vkjL = self.v[k][j][i].lb
                        vkjU = self.v[k][j][i].ub
                    else: # indica que a variável possui lb negativo
                        vkjL = - self.v[k][j][i][1].ub  # lb
                        vkjU = self.v[k][j][i][0].ub  # ub

                    self.model += (self.v[j][k][i][0] - self.v[j][k][i][1]) == self.x[j][i] - self.x[k][i] - self.instance.times[k][i], 'MCLinV({},{},{})'.format(j, k, i)
                    self.model += (self.v[k][j][i][0] - self.v[k][j][i][1]) == self.x[k][i] - self.x[j][i] - self.instance.times[j][i], 'MCLinV({},{},{})'.format(k, j, i)
                    self.model += (self.u[j][k][i][0] - self.u[j][k][i][1]) >= 0, 'MCLin1({},{},{})'.format(j, k, i)
                    self.model += (self.u[j][k][i][0] - self.u[j][k][i][1]) >= vjkL * self.y[k][j][i], 'MCLin2({},{},{})'.format(j, k, i)
                    self.model += (self.u[j][k][i][0] - self.u[j][k][i][1]) >= (self.v[j][k][i][0] - self.v[j][k][i][1]) + vjkU * self.y[k][j][i] - vjkU, 'MCLin3({},{},{})'.format(j, k, i)
                    self.model += (self.u[j][k][i][0] - self.u[j][k][i][1]) <= vjkU * self.y[k][j][i], 'MCLin4({},{},{})'.format(j, k, i)
                    self.model += (self.u[j][k][i][0] - self.u[j][k][i][1]) <= (self.v[j][k][i][0] - self.v[j][k][i][1]) + vjkL * self.y[k][j][i] - vjkL, 'MCLin5({},{},{})'.format(j, k, i)
                    self.model += (self.u[k][j][i][0] - self.u[k][j][i][1]) >= 0, 'MCLin6({},{},{})'.format(k, j, i)
                    self.model += (self.u[k][j][i][0] - self.u[k][j][i][1]) >= vkjL * self.y[j][k][i], 'MCLin7({},{},{})'.format(k, j, i)
                    self.model += (self.u[k][j][i][0] - self.u[k][j][i][1]) >= (self.v[k][j][i][0] - self.v[k][j][i][1]) + vkjU * self.y[j][k][i] - vkjU, 'MCLin8({},{},{})'.format(k, j, i)
                    self.model += (self.u[k][j][i][0] - self.u[k][j][i][1]) <= vkjU * self.y[j][k][i], 'MCLin9({},{},{})'.format(k, j, i)
                    self.model += (self.u[k][j][i][0] - self.u[k][j][i][1]) <= (self.v[k][j][i][0] - self.v[k][j][i][1]) + vkjL * self.y[j][k][i] - vkjL, 'MCLin10({},{},{})'.format(k, j, i)

        # constraints (5)
        for j in range(self.instance.n):
            self.model += self.c - self.x[j][self.instance.machines[j][self.instance.m - 1]] >= self.instance.times[j][
                self.instance.machines[j][self.instance.m - 1]], 'makespan({})'.format(j)
        self.model.write(
            '{}_modelMCLinNonNeg.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))


    # cutpool
    def optmizeCuts(self):
        # self.model.cuts_generator = Compact_CutPool(self.instance)
        self.model.optimize(max_seconds=1800)

    # relax model
    def relax(self):
        newConstraints = True
        self.model.relax()
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()

        while newConstraints:
            self.iterationsCuts += 1
            newConstraints = False
            self.model.verbose = 0
            self.model.optimize()
            self.clique_heuristic()
            if firstObjValue == 0:  # first execution
                firstObjValue = self.model.objective_value
                # execute clique (must be only executed at the beginning of the cuts)
                for a in range(self.instance.m):
                    # clique_cuts = 0
                    clique_cuts = self.clique_cuts_best(a)
                    self.totalCliqueCuts += clique_cuts

            # self.printSolution()
            # input('solution')
            hasCuts = 0
            # start = time()
            for a in range(self.instance.m):
                # print('Machine {}'.format(a))
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
                # print('For machine {} were found {} basic cuts, {} two jobs cuts, {} clique cuts, {} triangle cuts, '
                #       '{} basic cuts epsilon, {} half cuts and {} late jobs cuts'.format(a, basic_cuts, two_job_cuts,
                #                                                                          clique_cuts, triangle_cuts,
                #                                                                          basic_cuts_epsilon, half_cuts,
                #                                                                          late_jobs_cuts))
                hasCuts += basic_cuts + basic_cuts_epsilon + half_cuts + late_jobs_cuts + two_job_cuts + triangle_cuts
                # self.model.write('teste.lp')
                # input('iteracao completa')
            # end = time()
            # print('Time elapsed for iteration {}: {}s'.format(self.iterationsCuts, round(end - start, 2)))
            cutsFound += hasCuts
            # print('Cuts found: {}'.format(hasCuts))
            # print('objective value : {}'.format(self.model.objective_value))
            # for j in range(self.instance.n):
            #     for i in range(self.instance.m):
            #         print('x({},{}) = {}\t'.format(j, i, round(self.x[j][i].x, 2)), end='')
            #     print()
            # input()
            if hasCuts > 0:
                newConstraints = True
            gainObj = self.model.objective_value - firstObjValue
        # print(
        #     'Number of iterations: {}. Gain of objective value: {}. Total of cuts found: {}. Objective value: {}.'.format(
        #         self.iterationsCuts - 1, gainObj, cutsFound, self.model.objective_value))
        # print('Were found a total of {} basic cuts, {} two jobs cuts, {} clique cuts, {} triangle cuts, '
        #       '{} basic cuts epsilon, {} half cuts and {} late jobs cuts'.format(self.totalBasicCuts,
        #                                                                          self.totalTwoJobsCuts,
        #                                                                          self.totalCliqueCuts,
        #                                                                          self.totalTriangleCuts,
        #                                                                          self.totalBasicCutsEpsilon,
        #                                                                          self.totalHalfCuts,
        #                                                                          self.totalLateJobCuts))
        # self.model.write(
        #     '{}_relax_model.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))
        end = time()
        print('Time elapsed relaxation: {}s'.format(round(end - start, 2)))
        return

    # optimize model
    def optimizeInteger(self):
        self.model.optimize(max_seconds=1800)

    # optimize relax model
    def optimizeRelax(self):
        self.model.relax()
        self.model.optimize()

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
        m = self.m
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
        c = m.add_var(var_type=CONTINUOUS, lb=-100 * self.instance.K, name='C')
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')

        m.objective = c

        m += c == rs - ls, 'C'
        m += rs - ls <= 0, 'violado'
        m += rs == xsum(self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n)), 'RS'
        m += ls == xsum(self.instance.times[j][a] * e_p[j] for j in range(self.instance.n)) + xsum(
            self.instance.times[j][a] * self.instance.times[i][a] * xij[j][i] for i in range(self.instance.n) for j
            in range(i + 1, self.instance.n)), 'LS'

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
        if m.status != OptimizationStatus.OPTIMAL:
            return 0

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
        # m.write('teste_cut.lp')
        # print(c_name, m.objective_value)
        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     print(c_name, m.objective_value)
        #     print(c)
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        # self.model.remove(c)
        self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, c_name

        # self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, 'cut_basic_best' \
        #                                                                              '{}({},{})'.format(
        #     self.iterationsCuts, ''.join(str(i) for i in S), a)
        return 1

    def basic_cuts_plus_epsilon_best(self, a, k):
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        c = m.add_var(var_type=CONTINUOUS, lb=-100 * self.instance.K, name='C')
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')

        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                xij[i][j] = xij[j][i]

        var = xsum(
            self.instance.times[j][a] * self.instance.times[i][a] * xij[i][j] for i in range(self.instance.n) for j in
            range(i + 1, self.instance.n))  # part 2 of cut
        var += self.instance.e[k][a] * xsum(
            self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k)  # part 3 of cut
        var -= xsum(
            self.y[j][k][a].x * max(self.instance.e[k][a] - self.instance.e[j][a], 0) * self.instance.times[j][a] *
            x_aux[j] for j in range(self.instance.n) if
            j != k and not isinstance(self.y[j][k][a], int))  # part 4 of cut if i == j (just expand the multiplication)
        var -= xsum(
            self.y[i][k][a].x * max(self.instance.e[k][a] - self.instance.e[i][a], 0) * self.instance.times[j][a] *
            xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if
            i != k and j != k and i != j and not isinstance(self.y[i][k][a],
                                                            int))  # part 4 of cut if i != j (just expand the multiplication)

        m.objective = c
        m += c == ls - rs, 'C'
        m += c <= 0, 'violado'
        m += ls == xsum(
            self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in
            range(self.instance.n)), 'LS'  # part 1 of cut
        m += rs == var, 'RS'
        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += x_aux[k] == 0, 'not_in_S'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'

        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        c_name = 'cut_basic_plus_epsilon_best({},{},{})'.format(''.join(str(i) for i in S), a, k)
        # self.printSolution()
        # for j in range(self.instance.n):
        #     print('{}*{}*{} + '.format(self.instance.times[j][a], self.x[j][a].x, x_aux[j]), end='')
        # print()
        # print(c_name, m.objective_value)
        # m.write('teste_cut.lp')
        # input()

        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     print(c_name, m.objective_value)
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        # self.model.remove(c)

        # right side
        rs = self.instance.e[k][a] * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                rs += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
            self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * \
                    self.p(S, a) >= rs, c_name

        # self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
        #     self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,
        #                                                                                                  a) >= ls, 'cut_basic_plus_epsilon_best{}({},{},{})'.format(
        #     self.iterationsCuts, ''.join(str(i) for i in S), a, k)
        return 1

    def half_cuts_best(self, a, k):
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        t = [m.add_var(var_type=INTEGER, lb=0, name='t({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='o({})'.format(i)) for i in range(self.instance.n)]
        e = m.add_var(var_type=INTEGER, name='E')
        c = m.add_var(var_type=CONTINUOUS, lb=-100 * self.instance.K, name='C')
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')
        # C = m.add_var(name='C', lb=-100 * self.instance.K)

        m.objective = c

        m += ls == self.x[k][a].x, 'LS'
        m += rs == e + xsum(
            self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k and
            not isinstance(self.y[j][k][a], int)), 'RS'
        m += c == ls - rs, 'C'
        m += c <= 0, 'violado'
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
        # m += C + e + xsum(
        #     self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k and
        #         not isinstance(self.y[j][k][a], int)) == self.x[k][a].x

        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        if (m.objective_value) > -0.00000001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        c_name = 'cut_half_best({},{},{})'.format(k, ''.join(str(i) for i in S), a)
        # self.printSolution()
        # print(c_name, m.objective_value)
        # m.write('teste_cut.lp')
        # print(e.x, end='')
        # for j in range(self.instance.n):
        #     if j != k:
        #         print(' + {}*{}*{}'.format(self.y[j][k][a].x, self.instance.times[j][a], x_aux[j]), end='')
        # input()
        # c = self.model.constr_by_name(c_name)
        # # # print(c_name)
        # if c is not None:
        #     print(c)
        #     print(c_name, m.objective_value)
        # #     print('{} {}'.format(self.x[k][a].x, m.objective_value))
        # #     print(self.x[k][a].x - self.E(S, a) - xsum(self.y[j][k][a].x * self.instance.times[j][a] for j in S if j != k))
        #     self.printSolution()
        #     self.model.write('teste.lp')
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        #     # self.model.remove(c)
        self.model += self.x[k][a] - self.E(S, a) - xsum(
            self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= 0, c_name

        # self.model += self.x[k][a] - self.E(S, a) - xsum(
        #     self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= 0, 'cut_half_best{}({},{},{})'.format(
        #     self.iterationsCuts, k, ''.join(str(i) for i in S), a)
        # input()
        return 1

    def late_job_cuts_best(self, a, k, l):
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')
        C = m.add_var(var_type=CONTINUOUS, name='C', lb=-100 * self.instance.K)
        var = self.instance.e[l][a]
        var += xsum(self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n)
                    if j != k and not isinstance(self.y[j][k][a], int))
        var -= xsum(self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0) * x_aux[j]
                    for j in range(self.instance.n) if j != l and not isinstance(self.y[j][l][a], int))

        m.objective = C
        m += ls == self.x[k][a].x, 'LS'
        m += rs == var, 'RS'
        m += C == ls - rs, 'C'
        m += C <= 0, 'violado'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'

        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        c_name = 'cut_late_job_best({},{},{},{})'.format(''.join(str(i) for i in S), a, k, l)
        # print(c_name, m.objective_value)
        # self.printSolution()
        # m.write('teste_cut.lp')
        # input()
        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     print(c_name, m.objective_value)
        #     self.printSolution()
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        # self.model.remove(c)
        self.model += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
            self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
                    self.instance.e[l][a], c_name

        # self.model += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
        #     self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
        #             self.instance.e[l][a], 'cut_late_job_best{}({},{},{},{})'.format(self.iterationsCuts,
        #                                                                              ''.join(str(i) for i in S), a,
        #                                                                              k, l)
        return 1

    def two_jobs_cuts_best(self, a):
        # print('Best two jobs cuts')
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')
        c = m.add_var(var_type=CONTINUOUS, lb=-100 * self.instance.K, name='c')
        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                xij[i][j] = xij[j][i]

        var = 0
        var2 = 0
        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                var += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a].x * \
                       xij[i][j]
                var += (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a].x * \
                       xij[i][j]
                var2 += self.instance.times[i][a] * self.instance.times[j][a] * xij[i][j]
                var2 += self.instance.e[i][a] * self.instance.times[j][a] * xij[i][j]
                var2 += self.instance.e[j][a] * self.instance.times[i][a] * xij[i][j]
        m.objective = c

        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += xsum(x_aux[j] for j in range(self.instance.n)) == 2, 'must_be_2'
        m += ls == var, 'LS'
        m += rs == var2, 'RS'
        m += c == ls - rs, 'C'
        m += c <= 0, 'violado'

        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        if m.objective_value > -0.000001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        # apply cut to the model
        c_name = 'cut_two_jobs({},{},{})'.format(S[0], S[1], a)
        # m.write('teste_cut.lp')
        # print(c_name, m.objective_value)
        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     self.printSolution()
        #     i = S[0]
        #     j = S[1]
        #     print(c_name)
        #     print(c)
        #     print(m.objective_value)
        #     print('({} + {} - {}) * {} + ({} + {} - {}) * {} >= {}*{} + {}*{} + {}*{}'.format(
        #         self.instance.times[i][a], self.instance.e[i][a], self.instance.e[j][a], self.x[i][a].x,
        #         self.instance.times[j][a], self.instance.e[j][a], self.instance.e[i][a], self.x[j][a].x,
        #         self.instance.times[i][a], self.instance.times[j][a],
        #         self.instance.e[i][a], self.instance.times[j][a],
        #         self.instance.e[j][a], self.instance.times[i][a]
        #     ))
        #     m.write('teste_cut.lp')
            # self.model.write('teste.lp')
            # input()
            # self.model.remove(c)

        self.two_jobs_cuts(S[0], S[1], a)
        return 1

    def basic_cuts(self, S, a):
        cuts = 0

        # right side
        rs = 0
        for j in S:
            rs += self.instance.times[j][a] * self.x[j][a].x

        # left side

        ls = self.E(S, a) * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        # check violated cut
        if (ls - rs) > 0.00001:  # if rs < ls:
            cuts += 1
            self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, 'basic_cuts{}({},{})'.format(
                self.iterationsCuts,
                ''.join(str(i) for i in S), a)

        # reverse
        # right side
        # rs = 0
        # for j in S:
        #     rs += self.instance.times[j][a] * (self.c.x - self.x[j][a].x)

        # # left side
        # ls = self.F(S, a) * self.p(S, a)

        # for j in S:
        #     ls += self.instance.times[j][a] * self.instance.times[j][a]

        # for i in range(len(S)):
        #     for j in range(i + 1, len(S)):
        #         ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        # # check violated cut
        # if (ls - rs) > 0.00001: # if rs < ls:
        #     cuts += 1
        #     self.model += xsum(self.instance.times[j][a] * (self.c - self.x[j][a]) for j in
        #                        S) >= ls, 'basic_cuts_reverse{}({},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a)

        return cuts

    # two-job cuts
    def two_jobs_cuts(self, i, j, a):
        # right side
        rs = self.instance.times[i][a] * self.instance.times[j][a] + \
             self.instance.e[i][a] * self.instance.times[j][a] + \
             self.instance.e[j][a] * self.instance.times[i][a]
        # print('corte: {} {}'.format(rs, ls))

        self.model += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
            a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
                        a] >= rs, 'cut_two_jobs({},{},{})'.format(i, j, a)

        # self.model += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
        #     a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
        #                 a] >= ls, 'cut_two_jobs{}({},{},{})'.format(self.iterationsCuts, i, j, a)
        return 1

    def clique_heuristic(self, a, steps):
        x_bar = [self.x[j][a].x for j in range(self.instance.n)]
        p = [self.instance.times[j][a] for j in range(self.instance.n)]
        est = [self.instance.est[j][a] for j in range(self.instance.n)]
        clique = Clique(x_bar, p, est, steps)
        escolhidos, t, custos = clique.annealing()
        cuts = 0
        if escolhidos.ndim > 1:
            # print(escolhidos, t, custos)
            cuts = len(escolhidos)
            for e in range(len(escolhidos)):
                c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, e, a)
                self.model += xsum(t[e][i] * self.x[i][a] for i in range(self.instance.n)) >= 1, c_name
                # print('{} : Chosen = {} ; t = {} ; cost = {}'.format(c_name, escolhidos[e], t[e], custos[e]))
        elif escolhidos.ndim == 1 and escolhidos.size > 0:
            cuts = 1
            c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, 1, a)
            self.model += xsum(t[i] * self.x[i][a] for i in range(self.instance.n)) >= 1, c_name
            # print('{} : Chosen = {} ; t = {} ; cost = {}'.format(c_name, escolhidos, t, custos))
        # print(x_bar, p, est)
        # input('machine: {}'.format(a))
        # self.model.write('teste.lp')
        # input('lp escrito')
        return cuts

    def LAHC(self, a, steps, l):
        x_bar = [self.x[j][a].x for j in range(self.instance.n)]
        p = [self.instance.times[j][a] for j in range(self.instance.n)]
        est = [self.instance.est[j][a] for j in range(self.instance.n)]
        clique = Clique(x_bar, p, est, steps)
        best_chosen, best_t, best_custo = clique.LAHC(l)
        cuts = 0
        if best_custo < 100:
            cuts = 1
            c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, 1, a)
            self.model += xsum(best_t[0][i] * self.x[i][a] for i in range(self.instance.n)) >= 1, c_name
            # print('{} : Chosen = {} ; t = {} ; cost = {}'.format(c_name, escolhidos, t, custos))
        # print(x_bar, p, est)
        # input('machine: {}'.format(a))
        # self.model.write('teste.lp')
        # input('lp escrito')
        return cuts


    def clique_cuts(self, S, a):
        s_aux = list(range(len(S)))
        perms = np.asarray(list(permutations(s_aux)))
        perms = [np.asarray(s) for s in perms]

        K = [[0] * len(S) for i in range(len(perms))]

        for i in range(len(K)):
            soma = 0
            for j in range(len(K[i])):
                if j == 0:
                    soma = self.instance.est[S[perms[i][j]]][a]
                else:
                    soma = max(soma + self.instance.times[S[perms[i][j - 1]]][a], self.instance.est[S[perms[i][j]]][a])
                for aux in range(j, len(K[i])):
                    K[i][perms[i][aux]] = soma
        m = self.m
        m.clear()
        m.verbose = 0
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        # c = m.add_var(var_type=CONTINUOUS, lb=0, name='c')
        m.objective = xsum(self.x[S[i]][a].x * t[i] for i in range(len(S)))

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)
        # m += c == xsum(self.x[S[i]][a].x * t[i] for i in range(len(S))), 'C'
        # m += c <= 1, 'violado'
        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        soma = m.objective_value
        if soma > 0.99999999:
            return 0
        # x_bar = []
        # t_bar = []
        # for i in range(len(S)):
        #     x_bar.append(self.x[S[i]][a].x)
        #     t_bar.append(t[i].x)
        c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, ''.join(str(i) for i in S), a)

        # print(x_bar, t_bar, soma, c_name)
        # m.write('teste_cut.lp')
        # print(S)
        # input('clique machine {}'.format(a))

        self.model += xsum(round(t[i].x, 10) * self.x[S[i]][a] for i in range(len(S))) >= 1, c_name

        return 1


    """
    # sc = size of smaller clique
    # lc = size of larger clique
    # maxCliques = maximum cliques to be investigated
    """

    def clique_cuts_best(self, a, lc=None, hc=None):
        cliquesFound = 0
        start = time()
        timeLimit = 600  # 10 minutes
        lc = self.lc if lc is None else min(lc, self.instance.n)
        hc = self.hc if hc is None else min(hc, self.instance.n)
        for sizeS in range(hc, lc, -1):
            print('LC: {} HC: {} Machine: {}'.format(lc, sizeS, a))
            comb = combinations(list(range(0, self.instance.n)), sizeS)  # combinations of all possibles jobs of size sizeS
            comb = list(comb)
            dict = {}
            for s in range(len(comb)):
                end = time()
                if (end - start) >= timeLimit:
                    return 0
                dist = 0
                S = comb[s]
                for i in list(range(0, sizeS)):
                    for j in list(range(i + 1, sizeS)):
                        dist += abs(self.x[S[i]][a].x - self.x[S[j]][a].x)
                dict[s] = dist
            i = 0
            end = time()
            if (end - start) >= timeLimit:
                # print('Time limit for enumaration of cliques')
                return 0
            for key, value in sorted(dict.items(), key=lambda l: l[1]):
                end = time()
                if (end - start) >= timeLimit:
                    # print('Time limit for finding cliques. Cliques found: {}'.format(cliquesFound))
                    return cliquesFound
                S = comb[key]
                cliquesFound += self.clique_cuts(S, a)
                i += 1
                if i >= self.maxCliques:
                    # print('Max clique found. Cliques found: {}'.format(i))
                    break
        return cliquesFound

    # triangle cuts
    def triangle_cuts(self, i, j, k, a):
        soma = self.y[i][j][a].x + self.y[j][k][a].x + self.y[k][i][a].x
        cuts = 0
        if soma > 2:
            cuts = 1
            self.model += self.y[i][j][a] + self.y[j][k][a] + self.y[k][i][
                a] <= 2, 'triangle_cuts{}({},{},{},{})'.format(self.iterationsCuts, i, j, k, a)
        return cuts

    def triangle_cuts_best(self, a):
        # print('Best triangle cuts')
        m = self.m
        m.clear()
        m.verbose = 0
        xij = [[0 for i in range(self.instance.n)] for j in range(self.instance.n)]
        for i in range(self.instance.n):
            for j in range(self.instance.n):
                if i != j:
                    if not isinstance(self.y[i][j][a], int) or not isinstance(self.y[i][j][a], int):
                        xij[i][j] = m.add_var(var_type=BINARY, name='xij({},{})'.format(i, j))

        var = - xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if i != j
                     and isinstance(xij[i][j], Var))
        m.objective = var

        m += xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i
                  and isinstance(xij[i][j], Var)) >= 3, 'minimum_3'
        for i in range(self.instance.n):
            m += xsum(xij[i][j] for j in range(self.instance.n) if j != i and isinstance(xij[i][j], Var)) \
                 - xsum(xij[j][i] for j in range(self.instance.n) if j != i and isinstance(xij[i][j], Var)) == 0, \
                 'flow({})'.format(i)
        try:
            m += xsum(
                self.y[i][j][a].x * xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i
                and not isinstance(self.y[i][j][a], int) and isinstance(xij[i][j], Var)) \
                 - xsum(xij[i][j] for i in range(self.instance.n) for j in
                        range(self.instance.n) if j != i and isinstance(xij[i][j], Var)) >= -1, 'triangle_cut'
        except:
            print(self.model.status, self.node)
            for i in range(self.instance.n):
                for j in range(self.instance.n):
                    if i != j:
                        if isinstance(self.y[i][j][a], int):
                            print('{} {:.4f}'.format(self.y[i][j][a].name, self.y[i][j][a].x))
                        else:
                            print('y({},{},{}) ----'.format(i, j, a))
            input('erro')
        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        if m.objective_value > -0.000000001 or m.status == OptimizationStatus.INFEASIBLE:
            return 0

        S = []
        names = []
        soma = 0
        for i in range(self.instance.n):
            for j in range(self.instance.n):
                if not isinstance(xij[i][j], int):
                    if xij[i][j].x > 0.99999999 and not isinstance(self.y[i][j][a], int):
                        names.append('({}{})'.format(i, j))
                        S.append(self.y[i][j][a])
                        soma += self.y[i][j][a].x * xij[i][j].x

        if len(S) <= 1 or soma - (len(S) - 1) <= 0.000000001:
            return 0
        if len(names) < 3:
            m.write('teste_cut.lp')
            print(names, m.objective_value)
            input('deu merda')
        c_name = 'cut_triangle_best({},{})'.format(''.join(str(i) for i in names), a)
        # print(c_name)
        # m.write('teste_cut.lp')
        # input()
        # print(c_name)
        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     print(c_name)
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        # self.model.remove(c)

        self.model += xsum(i for i in S) <= len(S) - 1, c_name
        return 1

    def basic_cuts_plus_epsilon(self, S, a, k):
        # right side
        rs = 0
        for j in S:
            rs += self.instance.times[j][a] * self.x[j][a].x
        # left side
        ls = 0
        soma1 = 0
        for j in range(len(S)):
            for i in range(j + 1, len(S)):
                soma1 += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        soma1 += self.instance.e[k][a] * self.p(S, a)

        ls += soma1

        soma2 = 0
        for j in S:
            soma2 += self.y[j][k][a].x * max(self.instance.e[k][a] - self.instance.e[j][a], 0)
        ls -= soma2 * self.p(S, a)
        cuts = 0

        if (ls - rs) > 0.00001:  # if rs < ls:
            cuts = 1
            self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
                self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,
                                                                                                             a) >= soma1, 'basic_cuts_epsilon{}({},{},{})'.format(
                self.iterationsCuts, ''.join(str(i) for i in S), a, k)

        return cuts

    def half_cuts(self, S, a, k):
        # left side
        ls = self.E(S, a)

        for j in S:
            if j != k:
                ls += self.y[j][k][a].x * self.instance.times[j][a]

        cuts = 0
        if (ls - self.x[k][a].x) > 0.00001:  # self.x[k][a].x < ls:
            cuts = 1
            self.model += self.x[k][a] - xsum(
                self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= self.E(S,
                                                                                            a), 'half_cuts{}({},{},{})'.format(
                self.iterationsCuts, ''.join(str(i) for i in S), k, a)

        return cuts

    def late_job_cuts(self, S, a, k, l):
        # left side
        ls = self.instance.e[l][a]

        for j in S:
            if j != k:
                ls += self.y[j][k][a].x * self.instance.times[j][a]

        for j in S:
            if j != l:
                ls -= self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0)

        cuts = 0

        if (ls - self.x[k][a].x) > 0.00001:  # self.x[k][a].x < ls:
            cuts = 1
            self.model += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
                self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
                          self.instance.e[l][a], 'late_job_cuts{}({},{},{},{})'.format(self.iterationsCuts,
                                                                                       ''.join(str(i) for i in S), a, k,
                                                                                       l)
        return cuts

    def printSolution(self):
        print("C: ", self.c.x)
        for j in range(self.instance.n):
            for i in range(self.instance.m):
                print('x({},{}) = {:.2f}\t '.format(j, i, self.x[j][i].x), end='')
            print()


    def testCliqueMIP(self, lc=None, hc=None):
        newConstraints = True
        self.model.relax()
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()
        self.iterationsCuts = 0
        self.model.verbose = 0
        self.model.optimize()
        firstObjValue = self.model.objective_value
        while newConstraints:
            self.iterationsCuts += 1
            hasCuts = 0
            newConstraints = False
            for a in range(self.instance.m):
                clique_cuts = self.clique_cuts_best(a, lc, hc)
                hasCuts += clique_cuts
            cutsFound += hasCuts
            if hasCuts > 0:
                self.model.optimize()
                print('Added {} cuts'.format(hasCuts))
                self.model.write('teste.lp')
                newConstraints = True
            if cutsFound > 500 and self.iterationsCuts > 2:
                newConstraints = False

        end = time()
        lastObjValue = self.model.objective_value
        elapsedTime = round(end - start, 2)
        return firstObjValue, lastObjValue, elapsedTime, cutsFound

    def testCliqueSA(self, steps):
        newConstraints = True
        self.model.relax()
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()
        self.iterationsCuts = 0
        self.model.verbose = 0
        self.model.optimize()
        firstObjValue = self.model.objective_value

        while newConstraints:
            self.iterationsCuts += 1
            hasCuts = 0
            newConstraints = False
            self.model.verbose = 0
            for a in range(self.instance.m):
                clique_cuts = self.clique_heuristic(a, steps)
                print('Cliques on machine {} = {}'.format(a, clique_cuts))
                hasCuts += clique_cuts
            cutsFound += hasCuts

            if hasCuts > 0:
                self.model.optimize()
                print('Added {} cuts with obj = {}'.format(hasCuts, self.model.objective_value))
                self.model.write('teste.lp')
                newConstraints = True
            # if self.iterationsCuts > 15:
            #     newConstraints = False

        end = time()
        lastObjValue = self.model.objective_value
        elapsedTime = round(end - start, 2)
        return firstObjValue, lastObjValue, elapsedTime, cutsFound


    def testCliqueLAHC(self, steps, l):
        newConstraints = True
        self.model.relax()
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()
        self.iterationsCuts = 0
        self.model.verbose = 0
        self.model.optimize()
        firstObjValue = self.model.objective_value

        while newConstraints:
            self.iterationsCuts += 1
            hasCuts = 0
            newConstraints = False
            self.model.verbose = 0
            for a in range(self.instance.m):
                clique_cuts = self.LAHC(a, steps, l)
                print('Cliques on machine {} = {}'.format(a, clique_cuts))
                hasCuts += clique_cuts
            cutsFound += hasCuts

            if hasCuts > 0:
                self.model.optimize()
                print('Added {} cuts with obj = {}'.format(hasCuts, self.model.objective_value))
                self.model.write('teste.lp')
                newConstraints = True
            # if self.iterationsCuts > 15:
            #     newConstraints = False

        end = time()
        lastObjValue = self.model.objective_value
        elapsedTime = round(end - start, 2)
        return firstObjValue, lastObjValue, elapsedTime, cutsFound
