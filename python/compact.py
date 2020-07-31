# compact.py
import string
from builtins import list
from config import Config
from JSSPInstance import JSSPInstance

from mip.callbacks import CutPool
from mip.model import Model, xsum, maximize
from mip.entities import Var
from mip.constants import INTEGER, BINARY, CONTINUOUS, OptimizationStatus
from mip.constants import CutType

import sys
from time import process_time, time
import numpy as np
import copy

from compact_cutpool import SubTourCutGenerator
from SA_clique import Clique

from itertools import permutations, combinations, product


class Compact:
    def __init__(self, conf: Config):
        self.config = conf
        self.instance = JSSPInstance(conf.get_property("instance_name"), conf.get_property("horizon"))
        self.m = Model(solver_name="cbc")
        self.m.verbose = 0
        self.model = Model(name='compact', solver_name="cbc")
        self.model.verbose = 1
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
        self.cp = CutPool()


    # build the problem with big-M
    def constructProblemM(self):
        # self.instance.print()

        self.c = self.model.add_var(var_type=INTEGER, lb=0, ub=self.instance.K, name="C")
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                self.x[j][i] = self.model.add_var(var_type=INTEGER, lb=self.instance.est[j][i], ub=self.instance.lst[j][i], name='x(j{},m{})'.format(j, i))
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

    # build the problem with big-M
    def constructProblemMSubCycle(self):
        # self.instance.print()
        self.y = [
            [[0 for i in range(self.instance.m)] for k in range(self.instance.n+1)] for j in range(self.instance.n)]

        self.c = self.model.add_var(var_type=INTEGER, name="C")
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                self.x[j][i] = self.model.add_var(var_type=INTEGER, name='x(j{},m{})'.format(j, i))
                for k in range(j+1, self.instance.n):
                    self.y[j][k][i] = self.model.add_var(var_type=BINARY, name='y(j{},k{},m{})'.format(j, k, i))
                    self.y[k][j][i] = self.model.add_var(var_type=BINARY, name='y(j{},k{},m{})'.format(k, j, i))
                self.y[j][self.instance.n][i] = self.model.add_var(var_type=BINARY, name='y(j{},k{},m{})'.format(j, self.instance.n, i))

        self.model.objective = self.c

        # constraints (2)
        for j in range(self.instance.n):
            for i in range(1, self.instance.m):
                self.model += self.x[j][self.instance.machines[j][i]] - self.x[j][self.instance.machines[j][i - 1]] >= \
                              self.instance.times[j][self.instance.machines[j][i - 1]], 'ord({},{})'.format(j, i)

        # constraints (3-4)
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                for k in range(self.instance.n):
                    if k != j:
                        self.model += self.x[j][i] - self.x[k][i] + self.instance.K * (1 - self.y[j][k][i]) >= \
                                      self.instance.times[k][i], 'psy({},{},{})'.format(j, k, i)
                self.model += xsum(self.y[j][k][i] for k in range(self.instance.n + 1)) == 1, 'sum({},{})'.format(j, i)
            for k in range(self.instance.n+1):
                self.model += xsum(self.y[j][k][i] for j in range(self.instance.n) if j != k) <= 1, 'sum2({},{})'.format(k, i)



        # constraints (5)
        for j in range(self.instance.n):
            self.model += self.c - self.x[j][self.instance.machines[j][self.instance.m - 1]] >= self.instance.times[j][
                self.instance.machines[j][self.instance.m - 1]], 'makespan({})'.format(j)
        self.model.write(
            '{}_modelMSubCycle.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))

    def optimizeSubCycle(self):
        # self.model.cuts_generator = SubTourCutGenerator(self.instance, self.x, self.y)
        self.model.optimize()
        while self.subciclos() > 0:
            self.model.write('teste.lp')
            self.model.optimize()
        # self.printSolution()


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
                    print(self.v[j][k][i], self.v[k][j][i])
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
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()

        self.model.optimize(relax=True)
        self.model.verbose = 0
        self.totalBasicCuts = 0
        self.totalBasicCutsEpsilon = 0
        self.totalTriangleCuts = 0
        self.totalCliqueCuts = 0
        self.totalTwoJobsCuts = 0
        self.totalLateJobCuts = 0
        self.totalHalfCuts = 0

        firstObjValue = self.model.objective_value
        # self.printSolution()
        # if self.config.get_property('clique_cuts') == 1:
        #     for a in range(self.instance.m):
        #         # clique_cuts = 0
        #         clique_cuts = self.clique_cuts_best(a)
        #         self.totalCliqueCuts += clique_cuts
        self.iterationsCuts = 0
        while newConstraints:
            self.cp = CutPool()
            self.iterationsCuts += 1
            newConstraints = False
            # self.clique_heuristic()
            hasCuts = 0
            # start = time()
            print('iteration', self.iterationsCuts)
            # input()
            clique_cuts = 0
            if self.config.get_property('mip_general_cliques_sucessor') == 1:
                clique_cuts += self.executeSucessor()
                # self.totalCliqueCuts += clique_cuts

            if self.config.get_property('mip_general_cliques_not_sucessor') == 1:
                clique_cuts += self.executeNotSucessor()
                # self.totalCliqueCuts += clique_cuts

            if self.config.get_property('mip_general_cliques_quad') == 1:
                clique_cuts += self.executeQuad(3, 3)

            if self.config.get_property('mip_general_cliques_cuts') == 1:
                clique_cuts += self.mip_general_cliques()
                # self.totalCliqueCuts += clique_cuts
            print('cliques', clique_cuts)
            for a in range(self.instance.m):
                # print('machine',a)
                triangle_cuts = 0
                if self.config.get_property('triangle_cuts') == 1:
                    triangle_cuts = self.triangle_cuts_best(a)
                basic_cuts = 0
                if self.config.get_property('basic_cuts') == 1:
                    set_s = []
                    # basic_cuts += self.basic_cuts_best(a)
                    # basic_cuts += self.basic_cuts_reverse_best(a)
                    for i in range(2, self.instance.n + 1):
                        comb = list(combinations(list(range(self.instance.n)), i))
                        for S in comb:
                            cuts = self.basic_cuts(S, a)
                            basic_cuts += cuts
                            if cuts > 0:
                                set_s.append(S)
                                    
                if self.config.get_property('clique_cuts') == 1:
                    for s in set_s:
                        # clique_cuts = 0
                        if len(s) > 8:
                            continue
                        # print(s)
                        clique_cuts = self.clique_cuts(s, a)
                        self.totalCliqueCuts += clique_cuts
                
                two_job_cuts = 0
                if self.config.get_property('two_jobs_cuts') == 1:
                    two_job_cuts = self.two_jobs_cuts_best(a)

                late_jobs_cuts = 0
                half_cuts = 0
                basic_cuts_epsilon = 0
                for k in range(self.instance.n):
                    if self.config.get_property('half_cuts_best') == 1:
                        half_cuts += self.half_cuts_best(a, k)
                    if self.config.get_property('basic_cuts_plus_epsilon') == 1:
                        basic_cuts_epsilon = self.basic_cuts_plus_epsilon_best(a, k)
                    if self.config.get_property('late_job_cuts_cuts') == 1:
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
                hasCuts += basic_cuts + basic_cuts_epsilon + half_cuts + late_jobs_cuts + two_job_cuts + triangle_cuts + clique_cuts
                # self.model.write('teste.lp')
            # end = time()
            # print('Time elapsed for iteration {}: {}s'.format(self.iterationsCuts, round(end - start, 2)))
            cutsFound += hasCuts
            print('Cuts found: {}'.format(hasCuts))
            print('objective value : {}'.format(self.model.objective_value))
            # for j in range(self.instance.n):
            #     for i in range(self.instance.m):
            #         print('x({},{}) = {}\t'.format(j, i, round(self.x[j][i].x, 2)), end='')
            #     print()
            # input()
            if hasCuts > 0:
                newConstraints = True
            self.model.optimize(relax=True)
            self.model.write('teste.lp')
            # self.printSolution()


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
            # '{}_relax_model.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))
        # end = time()
        # print('Time elapsed relaxation: {}s'.format(round(end - start, 2)))
        return

    # optimize model
    def optimizeInteger(self):
        self.model.optimize(max_seconds=1800)

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
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j in range(self.instance.n)]
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
        m.write('teste_cut.lp')
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

    def basic_cuts_reverse_best(self, a):
        
        m = self.m
        m.clear()
        m.verbose = 1
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        f_p = [m.add_var(var_type=INTEGER, lb=0, name='f_p({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j in range(self.instance.n)]
        y = [m.add_var(var_type=INTEGER, lb=0, name='y({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='e({})'.format(i)) for i in range(self.instance.n)]
        z = m.add_var(var_type=INTEGER, name='E')
        c = m.add_var(var_type=CONTINUOUS, lb=-100 * self.instance.K, name='C')
        rs = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='rs')
        ls = m.add_var(var_type=INTEGER, lb=-100 * self.instance.K, name='ls')

        m.objective = c

        m += c == rs - ls, 'C'
        m += rs - ls <= 0, 'violado'
        m += rs == xsum(self.instance.times[j][a] * (self.c.x - self.x[j][a].x) * x_aux[j] for j in range(self.instance.n)), 'RS'
        m += ls == xsum(self.instance.times[j][a] * f_p[j] for j in range(self.instance.n)) + xsum(self.instance.times[j][a]*self.instance.times[j][a] * x_aux[j] 
            for j in range(self.instance.n)) + xsum(self.instance.times[j][a] * self.instance.times[i][a] * xij[j][i] for i in range(self.instance.n) for j
            in range(i + 1, self.instance.n)), 'LS'

        for j in range(self.instance.n):
            m += v[j] - self.instance.f[j][a] * x_aux[j] + self.instance.K * x_aux[
                j] == self.instance.K, 'eq26({})'.format(j)
            m += z - v[j] <= 0, 'eq27({})'.format(j)
            m += z - y[j] >= 0, 'eq28({})'.format(j)
            m += y[j] - self.instance.K * o[j] <= 0, 'eq29({})'.format(j)
            m += y[j] - v[j] <= 0, 'eq30({})'.format(j)
            m += y[j] - v[j] - self.instance.K * o[j] >= - self.instance.K, 'eq31({})'.format(j)
            # z*x(j)
            m += f_p[j] - self.instance.K * x_aux[j] <= 0, 'f_s1({})'.format(j)
            m += f_p[j] - z <= 0, 'f_s2({})'.format(j)
            m += f_p[j] - z - self.instance.K * x_aux[j] >= -self.instance.K, 'f_s3({})'.format(j)
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

        print(m.num_solutions, m.sol_pool_size)
        # input()
        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        # left side
        ls = self.F(S, a) * self.p(S, a)
        # print(self.F(S, a), '*', self.p(S, a), '+ ', end='')
        for i in range(len(S)):
            ls += self.instance.times[S[i]][a]*self.instance.times[S[i]][a]
            # print(self.instance.times[S[i]][a], '*', self.instance.times[S[i]][a], '+ ',  end='')
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
                # print(self.instance.times[S[j]][a], '*', self.instance.times[S[i]][a], '+ ', end='')

        c_name = 'cut_basic_reverse_best({},{})'.format(''.join(str(i) for i in S), a)
        m.write('teste_cut.lp')
        # print(c_name, m.objective_value)
        
        # c = self.model.constr_by_name(c_name)
        # if c is not None:
        #     print(c_name, m.objective_value)
        #     print(c)
        #     m.write('teste_cut.lp')
        #     input()
        #     return 0
        # self.model.remove(c)
        var = xsum(self.instance.times[j][a] * (self.c - self.x[j][a]) for j in S) >= ls
        # self.printSolution()
        print(var, c_name)
        
        # input()
        self.model += var, c_name

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
        ls = 0

        for j in S:
            ls += self.instance.times[j][a] * self.x[j][a].x

        # left side

        rs = self.E(S, a) * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                rs += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        
        # check violated cut
        if (rs - ls) > 0.00001:  # if rs < ls:
            cuts += 1
            # print(ls, rs, xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= rs)
            # input()
            self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= rs, 'basic_cuts{}({},{})'.format(
                self.iterationsCuts,
                ''.join(str(i) for i in S), a)

        # reverse
        # right side
        ls = 0
        for j in S:
            ls += self.instance.times[j][a] * (self.c.x - self.x[j][a].x)

        # left side
        rs = self.F(S, a) * self.p(S, a)

        for j in S:
            rs += self.instance.times[j][a] * self.instance.times[j][a]

        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                rs += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        # check violated cut
        
        if (rs - ls) > 0.00001: # if rs < ls:
            cuts += 1
            # print(ls, rs, xsum(self.instance.times[j][a] * (self.c - self.x[j][a]) for j in S) >= rs)
            linexpr = xsum(self.instance.times[j][a] * (self.c - self.x[j][a]) for j in S) >= rs
            if self.cp.add(linexpr):
                self.model += linexpr, 'basic_cuts_reverse{}({},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a)

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
        escolhidos, t, custos, minimum, maximum, exact = clique.annealing_mip()
        cuts = 0
        if t.ndim > 1:
            # print(escolhidos, t, custos)
            cuts = len(escolhidos)
            for e in range(len(escolhidos)):
                # print(''.join(str(int(i)) for i in escolhidos[e]))
                c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, ''.join(str(int(i)) for i in escolhidos[e]), a)
                self.model += xsum(t[e][i] * self.x[i][a] for i in range(self.instance.n)) >= 1, c_name
                # print('{} : Chosen = {} ; t = {} ; cost = {}'.format(c_name, escolhidos[e], t[e], custos[e]))
        elif t.ndim == 1 and t.size > 0:
            cuts = 1
            # print(escolhidos[0])
            # print(''.join(str(int(i)) for i in escolhidos[0]))
            c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, ''.join(str(int(i)) for i in escolhidos[0]), a)
            self.model += xsum(t[i] * self.x[i][a] for i in range(self.instance.n)) >= 1, c_name
            # print('{} : Chosen = {} ; t = {} ; cost = {}'.format(c_name, escolhidos, t, custos))
        # print(x_bar, p, est)
        # self.model.write('teste.lp')
        # input('lp escrito')
        return cuts, minimum, maximum, exact

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
        

        # for i in range(len(K)):
        #     soma = 0
        #     for j in range(len(K[i])):
        #         if j == 0:
        #             soma = self.instance.est[S[perms[i][j]]][a]
        #         else:
        #             soma = max(soma + self.instance.times[S[perms[i][j - 1]]][a], self.instance.est[S[perms[i][j]]][a])
        #         for aux in range(j, len(K[i])):
        #             K[i][perms[i][aux]] = soma
        for i in range(len(K)):
            pos = perms[i][0] 
            K[i][pos] = self.instance.est[S[pos]][a]
            for j in range(1,len(K[i])):
                pos = perms[i][j]
                pos_ant = perms[i][j-1]
                K[i][pos] = max(K[i][pos_ant] + self.instance.times[S[pos_ant]][a],  self.instance.est[S[pos]][a])
        m = self.m
        m.clear()
        m.verbose = 0
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        # c = m.add_var(var_type=CONTINUOUS, lb=0, name='c')
        m.objective = xsum(round(self.x[S[i]][a].x, 8) * t[i] for i in range(len(S)))
        # print('machine', a)
        # for i in S:
        #     print('i', ': est =', self.instance.est[i][a], 'p =', self.instance.times[i][a])
        # for i in range(len(K)):
        #     for j in range(len(K[i])):
        #         print('{}*{} '.format(K[i][j], t[j]), end='')
        #     print()
        # input()

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)
        # m += c == xsum(self.x[S[i]][a].x * t[i] for i in range(len(S))), 'C'
        # m += c <= 1, 'violado'
        m.optimize()

        if m.status != OptimizationStatus.OPTIMAL:
            return 0

        soma = round(sum(round(t[i].x, 10) * self.x[S[i]][a].x for i in range(len(S))), 8)
        if soma >= 0.99:
            return 0
        # x_bar = []
        # t_bar = []
        # for i in range(len(S)):
        #     x_bar.append(self.x[S[i]][a].x)
        #     t_bar.append(t[i].x)
        c_name = 'cut_clique_{}({},{})'.format(self.iterationsCuts, ''.join(str(i) for i in S), a)

        # print(x_bar, t_bar, soma, c_name)
        m.write('teste_cut.lp')
        # if (c_name == 'cut_clique_1(03478,5)'):
        #     for i in S:
        #         print(self.x[i][a], 'est = {} p = {}'.format(self.instance.est[i][a], self.instance.times[i][a]))
        #     for i in range(len(K)):
        #         for j in range(len(K[i])):
        #             print('{}*t({}) '.format(K[i][j], S[j]), end='')
        #         print()
        # print(S, soma)
        # print('obj:', xsum(round(self.x[S[i]][a].x, 8) * t[i] for i in range(len(S))))
        # m.write('teste_cut.lp')
        # print(c_name,':', xsum(round(t[i].x, 8) * self.x[S[i]][a] for i in range(len(S))))
        # input('clique machine {}'.format(a))

        self.model += xsum(round(t[i].x, 10) * self.x[S[i]][a] for i in range(len(S))) >= 0.99, c_name

        return 1


    """
    # sc = size of smaller clique
    # lc = size of larger clique
    # maxCliques = maximum cliques to be investigated
    """

    def clique_cuts_best(self, a, lc=None, hc=None):
        cliquesFound = 0
        # start = time()
        # timeLimit = 60000  # 10 minutes
        lc = self.lc if lc is None else min(lc, self.instance.n)
        hc = self.hc if hc is None else min(hc, self.instance.n)
        set_s = list()
        for sizeS in range(hc, lc-1, -1):
            # print('LC: {} HC: {} Machine: {}'.format(lc, sizeS, a))
            comb = combinations(list(range(0, self.instance.n)), sizeS)  # combinations of all possibles jobs of size sizeS
            comb = list(comb)
            set_s.extend(comb)
        
        # input(set_s)
        # for s in set_s:
        #     # end = time()
        #     # if (end - start) >= timeLimit:
        #     #     return 0
        #     dist = 0
        #     S = comb[s]
        #     for i in list(range(0, sizeS)):
        #         for j in list(range(i + 1, sizeS)):
        #             dist += abs(self.x[S[i]][a].x - self.x[S[j]][a].x)
        #     dict[s] = dist
        i = 0
        # end = time()
        # if (end - start) >= timeLimit:
        #     # print('Time limit for enumaration of cliques')
        #     return 0
        for S in set_s:
            # input(S)
            # end = time()
            # if (end - start) >= timeLimit:
            #     # print('Time limit for finding cliques. Cliques found: {}'.format(cliquesFound))
            #     return cliquesFound
            # S = comb[key]
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

        var = - xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if i != j and isinstance(xij[i][j], Var))
        m.objective = var

        m += xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i and isinstance(xij[i][j], Var)) >= 3, 'minimum_3'
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
            # input('erro')
        m.write('teste_cut.lp')
        # input('teste')
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
            # input('deu merda')
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
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    print('y({},{},{}) = {:.2f}\t '.format(j, k, i, self.y[j][k][i].x), end='')
                    print('y({},{},{}) = {:.2f}\t '.format(k, j, i, self.y[k][j][i].x), end='')
                    print()

    def testCliqueMIP(self, lc=None, hc=None):
        newConstraints = True
        # self.model.relax()
        gainObj = 0
        cutsFound = 0
        firstObjValue = 0
        start = time()
        self.iterationsCuts = 0
        self.model.verbose = 0
        self.model.optimize(relax=True)
        # self.printSolution()
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
                self.model.optimize(relax=True)
                # self.printSolution()
                print('Added {} cuts'.format(hasCuts))
                self.model.write('teste.lp')
                newConstraints = True
            # if cutsFound > 500 and self.iterationsCuts > 2:
            #     newConstraints = False

        end = time()
        lastObjValue = self.model.objective_value
        elapsedTime = round(end - start, 2)

        self.model.write('teste.lp')
        print(firstObjValue, lastObjValue, elapsedTime, cutsFound)
        # input()
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
        minimum = 0
        maximum = 0
        exact = 0
        maxTime = 60*60*3   # 3 hours
        initial = time()

        while newConstraints:
            self.iterationsCuts += 1
            hasCuts = 0
            newConstraints = False
            self.model.verbose = 0
            for a in range(self.instance.m):
                clique_cuts, minimum, maximum, exact = self.clique_heuristic(a, steps)
                end = time()
                if (end - initial) > maxTime:
                    break
                print('Cliques on machine {} = {}. Minimum = {}, maximum = {}, exacts found = {}. Elapsed: {:>4.3g}s'.format(a, clique_cuts, minimum, maximum, exact, end-initial))
                hasCuts += clique_cuts
            cutsFound += hasCuts

            end = time()
            if (end - initial) > maxTime:
                break

            if hasCuts > 0:
                self.model.optimize()
                print('Added {} cuts with obj = {}'.format(hasCuts, self.model.objective_value))
                # self.model.write('teste.lp')
                newConstraints = True

            if (end - initial) > maxTime:
                break

            # self.printSolution()
            # input()
            # if self.iterationsCuts > 15:
            #     newConstraints = False

        end = time()
        lastObjValue = self.model.objective_value
        elapsedTime = round(end - start, 2)
        return firstObjValue, lastObjValue, elapsedTime, cutsFound, minimum, maximum, exact


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

    def subciclos(self):
        yf = self.model.translate(self.y)
        V_ = range(self.instance.n)

        cp = CutPool()
        for a in range(self.instance.m):
            edges = [[] for i in range(self.instance.n)]
            for (u, v) in [(k, l) for (k, l) in product(V_, range(self.instance.n+1)) if k != l and abs(yf[k][l][a].x) > 1e-8]:
                print(yf[u][v][a].name, yf[u][v][a].x)
                edges[u].append(v)
            # print(edges)
            # input()
            cycles = []
            for i in range(self.instance.n):
                c = self.cycles(edges[i], i, edges, [i], [])
                if c:
                    for cycle in c:
                        cycles.append(cycle)
            # print('cycles', cycles, 'a ', a)
            # input()
            for c in cycles:
                soma = 0
                rhs = len(c) - 1
                for i in range(rhs):
                    u = c[i]
                    v = c[i+1]
                    soma += yf[u][v][a].x
                if soma > rhs - 1 + 1e-8:
                    cut = xsum(self.y[c[i]][c[i+1]][a] for i in range(rhs)) <= rhs - 1
                    # print(cut)
                    # input()
                    cp.add(cut)
        for cut in cp.cuts:
            print(cut)
            self.model += cut
        print('Total cuts: {}'.format(len(cp.cuts)))
        # input()
        return len(cp.cuts)

    def cycles(self, neighbors, first, edges, cycle, total):
        # print('neighbors', neighbors, 'first', first, 'edges', edges, 'cycle', cycle, 'cycles', total)
        if neighbors:
            # print(neighbors)
            for n in neighbors:
                if n >= len(edges) or n < first:  # if is the last one or already checked some
                    return
                # print(cycle)
                # print('edges[{}]'.format(n), edges[n], 'first', first, 'edges', edges, 'cycle', cycle, 'cycles', total)
                # input()
                if n == first:
                    c = cycle.copy()
                    c.append(n)
                    total.append(c)
                elif n not in cycle:
                    cycle.append(n)
                    self.cycles(edges[n], first, edges, cycle, total)
                    cycle.pop()
        # print('total', total)
        # input()
        return total

    def mip_generate_cut_class(self):
        cortes = 0
        if self.config.get_property("cut_type_gomory") == 1:
            cp = self.model.generate_cuts([CutType.GOMORY], 8192, 1e-4)
            for c in cp.cuts:
                self.iterationsCuts += 1
                if len(c.expr) < self.config.get_property("cut_type_gomory_limit"):
                    self.model.add_constr(c, 'gomory({})'.format(self.iterationsCuts))
            cortes += len(cp.cuts)
        if self.config.get_property("cut_type_zero_half") == 1:
            cp = self.model.generate_cuts([CutType.ZERO_HALF], 8192, 1e-4)
            for c in cp.cuts:
                self.iterationsCuts += 1
                self.model.add_constr(c, 'zero_half({})'.format(self.iterationsCuts))
            cortes += len(cp.cuts)
        if self.config.get_property("cut_type_mir") == 1:
            cp = self.model.generate_cuts([CutType.MIR], 8192, 1e-4)
            for c in cp.cuts:
                self.iterationsCuts += 1
                self.model.add_constr(c, 'mir({})'.format(self.iterationsCuts))
            cortes += len(cp.cuts)
        self.model.write('teste.lp')
        return cortes
        # self.model.optimize(relax=True)
        # self.printSolution()
        # input()


    def mip_general_cliques(self):
        # self.model.relax()
        self.model.verbose = 0
        self.model.optimize(relax=True)
        # self.printSolution()
        # self.model.write('teste.lp')
        # input()
        select = self.config.get_property("mip_general_cliques_select")
        param = self.config.get_property("mip_general_cliques_parameter")
        d = {}
        if select == 0:  # select_y
            d = self.select_tuples_y(param)
        elif select == 1:  # select_intersec
            d = self.select_tuples_intersec(param)
        elif select == 2:  # select_l
            d = self.select_tuples_l(param)
        elif select == 3:
            d = self.select_tuples_l_intersec(param)
        # print(d)
        # input('check')
        cortes = 0
        cliques = 0
        mip = 0
        if len(d) > 0:
            if select == 3:
                index = 0
                while cliques == 0 and index < 5:
                    cliques = self.solveGeneralClique(d[index])
                    index += 1
                    # print(cliques, index)
                    # input()
            else:
                cliques = self.solveGeneralClique(d)
            # cliques = self.general_cliques(d)
        total = cliques

        # mip = self.mip_generate_cut_class()
        mip = 0
        cortes = cliques + mip
        qtd = 0
        mip_maximum = self.config.get_property("cut_mip_maximum")
        self.iterationsCuts += 1
        if mip_maximum is None:
            mip_maximum = 0
        while cortes > 0:
            qtd += 1
            # self.model.relax()
            self.model.optimize(relax=True)
            # self.printSolution()
            # input()
            # self.printSolution()
            # self.model.write('teste.lp')
            # print('cortes: ', cortes)
            cliques = 0
            mip = 0
            # input()
            # input('teste.lp')
            if select == 0:  # select_y
                d = self.select_tuples_y(param)
            elif select == 1:  # select_intersec
                d = self.select_tuples_intersec(param)
            elif select == 2:  # select_l
                d = self.select_tuples_l(param)
            elif select == 3:  # select_l_intersec
                d = self.select_tuples_l_intersec(param)
            if len(d) > 0:
                if select == 3:
                    index = 0
                    cliques = 0
                    while cliques == 0 and index < 5:
                        cliques = self.solveGeneralClique(d[index])
                        index += 1
                        # print(cliques, index)
                        # input()
                else:
                    cliques = self.solveGeneralClique(d)
                # cliques = self.general_cliques(d)
                total += cliques
            # print(d)
            # input('check')

            # if self.iterationsCuts < mip_maximum:
            #     mip = self.mip_generate_cut_class()
            cortes = cliques + mip
            # print("cliques: ", cliques, 'mip:', mip)
        self.iterationsCuts = qtd
        print('iterações:', qtd, 'cortes cliques encontrados:', self.totalCliqueCuts, 'mip:', self.iterationsCuts, 'objective: ', self.model.objective_value)
        
        return qtd
        # input()
        # self.model.optimize()
        # self.model.write('teste.lp')
        # self.printSolution()
        # input()

    def select_tuples_random(self, max):
        if max > self.instance.n * self.instance.m:
            max = self.instance.n * self.instance.m
        d = []
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                d.append((j, i))
        test = []
        for i in range(max):
            r = np.random.randint(0, len(d))
            test.append(d[r])
            d.pop(r)
        return test

    def select_tuples_intersec(self, max):
        d = set()
        for i in range(self.instance.m):
            first = False
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    if self.x[j][i].x < self.x[k][i].x:
                        if (self.x[k][i].x - self.x[j][i].x) < self.instance.times[j][i]:
                            # print(self.x[k][i].x, self.x[j][i].x, self.instance.times[j][i])
                            if first == False:
                                first = True
                                if len(d) < max:
                                    d.add((j, i))
                                    # print(len(d))
                                    # input((j, i))
                                else:
                                    break
                            if len(d) < max:
                                d.add((k, i))
                                # print(len(d))
                                # input((k, i))
                            else:
                                break
                    else:
                        if (self.x[j][i].x - self.x[k][i].x) < self.instance.times[k][i]:
                            print(self.x[j][i].x, self.x[k][i].x, self.instance.times[k][i])
                            if first == False:
                                first = True
                                if len(d) < max:
                                    d.add((j, i))
                                    # print(len(d))
                                    # input((j, i))
                                else:
                                    break
                            if len(d) < max:
                                d.add((k, i))
                                # print(len(d))
                                # input((k, i))
                            else:
                                break
        return list(d)

    def select_tuples_y(self, max):
        d = set()
        set_y = {}
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    if self.y[j][k][i].x < 1 - 1e-8:
                        set_y[(j, k, i)] = abs(self.y[j][k][i].x - 0.5)
        set_y = sorted(set_y.items(), key=lambda l: l[1])
        i = 0
        
        # input()
        if max > len(set_y):
            max = len(set_y)

        # print(max)
        while i < max:
            key, value = set_y[i]
            # print(key, value)
            if (key[0], key[2]) not in d:
                i += 1
                d.add((key[0], key[2]))
            if i < max:
                if (key[1], key[2]) not in d:
                    i += 1
                    d.add((key[1], key[2]))
        
        return list(d)

    def select_tuples_l(self, l):
        d = set()
        set_aux = {}
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    if self.y[j][k][i].x < 1 - 1e-8:
                        set_aux[(j, k, i)] = abs(self.y[j][k][i].x - 0.5)
        set_aux = sorted(set_aux.items(), key=lambda l: l[1])
        i = 0
        # print(set_aux, len(set_aux))
        # input()

        key, value = set_aux[0]
        # print(key, value)

        # add two sides of most fractionary y
        d.add((key[0], key[2]))
        d.add((key[1], key[2]))

        machines = set()
        machines.add(key[2])

        # add the machines before (j, i) and (k, i) if exists
        j = key[0]
        k = key[1]
        i = key[2]
        order_j = self.instance.o[j][i]
        order_k = self.instance.o[k][i]
        # print(order_j, order_k)
        if order_j > 0:
            mach = self.instance.machines[j][order_j - 1]
            # print('j', j, 'mach', mach)
            machines.add(mach)
            d.add((j, mach))
        if order_k > 0:
            mach = self.instance.machines[k][order_k - 1]
            machines.add(mach)
            # print('k', k, 'mach', mach)
            d.add((k, mach))

        # print('d', d)
        # print('machines', machines)
        # add jobs j that aren't in d which machines i are in set machines that have the lowest value x_ji
        del set_aux
        set_aux = {}
        for i in machines:
            for j in range(self.instance.n):
                set_aux[(j, i)] = self.x[j][i].x
        set_aux = sorted(set_aux.items(), key=lambda l: l[1])
        # print(set_aux)

        i = 0
        max = l + len(d)
        while len(d) < max and i < len(set_aux):
            key, value = set_aux[i]
            # print(key, value)
            d.add((key[0], key[1]))
            i = i + 1
        # input(d)

        return list(d)

    def select_tuples_l_intersec(self, l):
        set_d = [list() for _ in range(5)]
        
        set_aux_y = {}
        for i in range(self.instance.m):
            for j in range(self.instance.n):
                for k in range(j + 1, self.instance.n):
                    if self.y[j][k][i].x < 1 - 1e-8:
                        set_aux_y[(j, k, i)] = abs(self.y[j][k][i].x - 0.5)
        set_aux_y = sorted(set_aux_y.items(), key=lambda l: l[1])
        i = 0
        for aux in range(5):
            d = set()
            key, value = set_aux_y[aux]
            # print(key, value)

            # add two sides of most fractionary y
            d.add((key[0], key[2]))
            d.add((key[1], key[2]))

            machines = set()
            machines.add(key[2])

            # add the machines before (j, i) and (k, i) if exists
            j = key[0]
            k = key[1]
            i = key[2]
            order_j = self.instance.o[j][i]
            order_k = self.instance.o[k][i]
            # print(order_j, order_k)
            if order_j > 0:
                mach = self.instance.machines[j][order_j - 1]
                # print('j', j, 'mach', mach)
                machines.add(mach)
                d.add((j, mach))
            if order_k > 0:
                mach = self.instance.machines[k][order_k - 1]
                machines.add(mach)
                # print('k', k, 'mach', mach)
                d.add((k, mach))

            # print('d', d)
            # print('machines', machines)
            # add jobs j that aren't in d which machines i are in set machines that have the lowest value x_ji
            set_aux = {}
            # self.printSolution()
            for i in machines:
                for j in range(self.instance.n):
                    for k in range(j+1, self.instance.n):
                        a0 = self.x[j][i].x
                        a1 = self.x[j][i].x + self.instance.times[j][i]
                        b0 = self.x[k][i].x
                        b1 = self.x[k][i].x + self.instance.times[k][i]
                        if (b0 > a1) or (a0 > b1):
                            continue
                        else:
                            o0 = max(a0, b0)
                            o1 = min(a1, b1)
                            if (o1 - o0) < 1e-8:
                                continue
                            set_aux[(j, k, i)] = o1 - o0
            set_aux = sorted(set_aux.items(), key=lambda lam: lam[1], reverse=True)
            # input(set_aux)
            index = 0
            maximum = l + len(d)
            while len(d) < maximum and index < len(set_aux):
                key, value = set_aux[index]
                index += 1
                j = key[0]
                k = key[1]
                i = key[2]
                d.add((j, i))
                if len(d) < l:
                    d.add((k, i))
            set_d[aux] = list(d)
        print(set_d)
        # input()
        return set_d

    def createK(self, S, pos, sol, K, last_job, last_machine, lenK, maxK):
        # print('S', S, 'len', len(S))
        # print('sol', sol)
        # print('pos', pos)
        # print('K', K[:lenK])
        # print('last_job', last_job)
        # print('last_machine', last_machine)
        # print('lenK', lenK)
        # print('maxK', maxK)
        # input()
        if len(S) < 1:  # caminho final válido
            # print(sol)
            # input('fim do caminho')

            if not any((K[:lenK]==sol).all(1)):
                # print('teste')
                K[lenK] = sol
                lenK = lenK + 1
                if lenK == len(K):
                    K = np.vstack((K, np.zeros(shape=(maxK, len(S)))))
            # input()
            return lenK
        
        for k in range(len(S)):
            # print('removing pos', k, 'from', S)
            (j, i) = S.pop(k)
            est_last_machine = 0
            est_last_job = 0
            # print((j, i))
            last_m = last_machine[j]
            # print('last_m', last_m)
            order_now = -1
            order_last = -1
            if last_m >= 0: # check for order conflict if a machine has been added for the job j
                order_now = self.instance.o[j][i]
                order_last = self.instance.o[j][last_m]
                # print(order_now, order_last, order_now < order_last)
                est_last_machine = sol[pos[(j, last_m)]] + self.instance.distances[j][last_m][i]

            if order_now < order_last:  # conflict found. Reinsert and go back
                # print('error. reinsert pos', k, 'from', S, 'and go back')
                S.insert(k, (j, i))
                # print('S', S)
                return lenK
            else: 
                last_j = last_job[i]
                # print('last_j', last_j)
                if last_j >= 0:
                    est_last_job = sol[pos[(last_j, i)]] + self.instance.times[last_j][i]
                est = max(self.instance.est[j][i], est_last_machine, est_last_job)
                # print(self.instance.est[j][i], est_last_machine, est_last_job, 'est', est)
                index = pos[(j, i)]
                est_antigo = sol[index]
                sol[index] = est
                last_machine[j] = i
                last_job[i] = j
                lenK = self.createK(S, pos, sol, K, last_job, last_machine, lenK, maxK)
                # print('reinsert pos', k, 'from', S, 'and go back')
                S.insert(k, (j, i))
                sol[index] = est_antigo
                last_machine[j] = last_m
                last_job[i] = last_j
                # print('S', S)
                # print('sol', sol)
                # input()
        return lenK

    def solveGeneralClique(self, S, name='general'):
        # print(S)
        P = list(permutations(range(len(S))))
        maxK = 100000
        lenK = 0
        pos = {S[i]: i for i in range(len(S))}
        K = np.zeros(shape=(maxK, len(S)))
        path = np.full(len(S), -1)
        last_job = np.full(self.instance.m, -1)  # last job added for machine i
        last_machine = np.full(self.instance.n, -1)  # last machine added in job j

        lenK = self.createK(S, pos, path, K, last_job, last_machine, lenK, maxK)
        # input(K[:lenK])

        m = self.m
        m.clear()
        m.verbose = 0
        t = {}
        for (j, i) in S:
            # print((j, i))
            t[(j, i)] = m.add_var(name='t({},{})'.format(j, i), lb=0, var_type=CONTINUOUS)

        # for a in t:
        #     print(a, t[a])
        qtd = 0
        # input()
        for p in K[:lenK]:
            # print(p)
            m += xsum(p[pos[(j, i)]]*t[(j,i)] for (j, i) in S) >= 1, 'c({})'.format(qtd)

            # for (j, i), est in p.items():
            #     # print(est[(j, i)],  t[(j, i)], end='')
            #     var += est * t[(j, i)]
            # # print(var)
            # m += var >= 1, 'c({})'.format(qtd)
            qtd += 1
        # input(K)
        # for (j,i) in d:
        #     print(self.x[j][i].x, t[(j, i)])
        # input()
        m.objective = xsum(self.x[j][i].x * t[(j, i)] for (j, i) in S)
        # m.write('teste_cuts.lp')
        # input('feito')
        m.optimize()
        if m.status != OptimizationStatus.OPTIMAL or m.objective_value > (1 - 1e-8):
            # if m.status == OptimizationStatus.OPTIMAL:
            #     print('obj: ', m.objective_value)
            # print('erro')
            return 0
        # print('obj: ', m.objective_value)
        # input()

        # print('solutions: ', m.num_solutions)
        # input()
        c_name = '{}_clique_cut_{}'.format(name, self.totalCliqueCuts)
        self.totalCliqueCuts += 1

        var = 0
        values = 'b;{}'.format(m.objective_value)
        for (job, machine) in S:
            var += self.x[job][machine] * round(t[(job, machine)].x, 12)
            if t[(job, machine)].x > 0.000001:
                values += ';{}; {}; {}; {}'.format(self.x[job][machine].name, self.x[job][machine].x, t[(job, machine)].name, t[(job, machine)].x)
        expr = var >= 1
        
        if self.cp.add(expr):
            self.model += expr, c_name
            print('{}: {}'.format(c_name, expr))
            print('{}: {}'.format(c_name, values))
            return 1
        return 0


    def general_cliques(self, d):
        # print('general_cliques')
        # print(d)
        pool = self.possible_paths(len(d), d, ([], {}), [])
        pool = [dict(s) for s in set(frozenset(d.items()) for d in pool)]
        # print(pool)
        # for p in pool:
        #     print(p)
        # input()

        m = self.m
        m.clear()
        m.verbose = 0
        t = {}
        for (j, i) in d:
            # print((j, i))
            t[(j, i)] = m.add_var(name='t({},{})'.format(j, i), lb=0, var_type=CONTINUOUS)

        # for a in t:
        #     print(a, t[a])
        qtd = 0
        for p in pool:
            # print(p)
            var = 0
            for (j, i), est in p.items():
                # print(est[(j, i)],  t[(j, i)], end='')
                var += est * t[(j, i)]
            # print(var)
            m += var >= 1, 'c({})'.format(qtd)
            qtd += 1
        # input(K)
        # for (j,i) in d:
        #     print(self.x[j][i].x, t[(j, i)])
        # input()
        m.objective = xsum(self.x[j][i].x * t[(j, i)] for (j, i) in d)
        m.write('teste_cuts.lp')
        # input()
        m.optimize()
        if m.status != OptimizationStatus.OPTIMAL or m.objective_value > (1 - 1e-8):
            if m.status == OptimizationStatus.OPTIMAL:
                print('obj: ', m.objective_value)
            print('erro')
            return 0
        print('obj: ', m.objective_value)
        # input()

        print('solutions: ', m.num_solutions)
        # input()
        c_name = 'cut_clique_{}'.format(self.totalCliqueCuts)
        self.totalCliqueCuts += 1

        var = 0
        for (job, machine) in d:
            var += self.x[job][machine] * t[(job, machine)].x
        print(var)
        self.model += var >= 1, c_name
        return 1


    def possible_paths(self, size, d, sol, pool):
        if len(sol[0]) == size:  # caminho final válido
            # est = self.calculate_est(sol)
            pool.append(sol[1].copy())
            # print('sol =', sol)
            # print('pool: ', pool)
            # input()
            return pool

        for k in range(len(d)):
            (j, i) = d.pop(k)
            if self.checkpath(sol, j, i):  # caso caminho seja válido até o moemnto
                # print((j, i), sol)
                # input()
                # sol[len(sol):] = [(j, i)]  # inserir na última posição
                pool = self.possible_paths(size, d, sol, pool)
                # print(d, sol)
                # input()
                sol[0].pop()  # remover da solução para colocar outra
                del sol[1][(j, i)]
                d.insert(k, (j, i))  # voltar a tupla para a posição original para não atrapalhar a iteração de d
            else:
                # print('faltam =', d, 'sol =', sol, 'rejeitado = (', j, ',', i, ')')
                # input()
                d.insert(k, (j, i))  # voltar a tupla para a posição original para não atrapalhar a iteração de d
                return pool # retornar pool (não precisa continuar a olhar o resto)
            # print('faltam =', d, 'sol =', sol)
            # input()

        return pool

    def checkpath(self, solution, j, i):
        est = self.instance.est[j][i]
        found_job = False
        found_machine = False
        # input(solution)
        for k in range(len(solution[0])-1, -1, -1):
            # print(k)
            order = self.instance.o[j][i]
            j_aux, i_aux = solution[0][k][0], solution[0][k][1]
            if i_aux == i and not found_machine:
                # print('found_machine', 'comparado = ', (j, i), 'no conjunto = ', (j_aux, i_aux), est, solution[1][(j_aux, i_aux)] + self.instance.times[j_aux][i_aux])
                est = max(est, solution[1][(j_aux, i_aux)] + self.instance.times[j_aux][i_aux])
                found_machine = True
            if j_aux == j:
                o = self.instance.o[j_aux][i_aux]
                if o > order:
                    # print('teste: ', (j, i), 'ordem: ', order, ' falhou em: ', solution[0][k], 'ordem: ', o)
                    return False
                # print('found_job', 'comparado = ', (j, i), 'no conjunto = ', (j_aux, i_aux), est, solution[1][(j_aux, i_aux)] + self.instance.times[j_aux][i_aux])
                est = max(est, solution[1][(j_aux, i_aux)] + self.instance.distances[j][i][i_aux])
                found_job = True
            if found_machine and found_job:
                break
        # print('est: ', est, 'lst', self.instance.lst[j][i])
        if est > self.instance.lst[j][i]:
            # print('est invalido: ', est, self.instance.lst[j][i])
            return False
        solution[0][len(solution[0]):] = [(j, i)]
        solution[1][(j, i)] = est
        return True

    def calculate_est(self, solution):
        est = {(j, i): self.instance.est[j][i] for (j, i) in solution}
        # input(est)
        for k in range(len(solution)):
            i = solution[k][1]
            j = solution[k][0]
            found_job = False
            found_machine = False
            for h in range(k+1, len(solution)):
                a = solution[h][1]
                b = solution[h][0]
                if not found_machine and a == i:
                    # print(solution[k], solution[h], est[solution[k]] + self.instance.times[j][i], est[solution[h]])
                    est[solution[h]] = max(est[solution[k]] + self.instance.times[j][i], est[solution[h]])
                    found_machine = True
                if not found_job and b == j:
                    # print('distance ', i, '->', a, '= ', self.instance.distances[j][i][a])
                    # print(solution[k], solution[h], est[solution[k]] + self.instance.distances[j][i][a], est[solution[h]])
                    est[solution[h]] = max(est[solution[k]] + self.instance.distances[j][i][a], est[solution[h]])
                    found_job = True
                if found_machine and found_job:
                    break
        # print('final est: ', est)
        # input()
        return est

    def executeSucessor(self):
        cliques = 0
        for i in range(self.instance.m - 1):
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    S = [(j, i), (k, i)]
                    o = self.instance.o[j][i]
                    if o < self.instance.m-1:
                        ii = self.instance.machines[j][o+1]
                        S.append((j, ii))
                    o = self.instance.o[k][i]
                    if o < self.instance.m-1:
                        ii = self.instance.machines[k][o+1]
                        S.append((k, ii))
                    cliques += self.solveGeneralClique(S, 'sucessor')
        return cliques

    def executeNotSucessor(self):
        cliques = 0
        for i in range(self.instance.m - 1):
            for j in range(self.instance.n):
                for k in range(j+1, self.instance.n):
                    for h in range(i+1, self.instance.m):
                        o1 = self.instance.o[j][i]
                        o2 = self.instance.o[j][h]
                        o3 = self.instance.o[k][i]
                        o4 = self.instance.o[k][h]
                        if abs(o1 - o2) > 1 and abs(o3 - o4) > 1:
                            S = [(j, i), (k, i), (j, h), (k, h)]
                            cliques += self.solveGeneralClique(S, 'notsucessor')
        return cliques

    # execute quadrant separation
    def executeQuad(self, sizeQ, sizeJ):
        cliques = 0
        for i in range(self.instance.m-sizeQ+1):
            for j in range(self.instance.n-sizeJ+1):
                S = []
                for k in range(j, j+sizeJ):
                    for h in range(i, i+sizeQ):
                        m2 = self.instance.machines[k][h]
                        S.append((k, m2))
                cliques += self.solveGeneralClique(S, 'quad')
        # self.model.write('teste.lp')
        # for c in self.cp.cuts:
        #     print(c)
        # input('acabou')
        return cliques

    def splitCuts(self):
        grid = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        delta = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
        continueSeparation = True
        
        flag = 0
        iterations = 0
        cuts = 0
        self.model.optimize(relax=True)
        print('solucao relaxada: {}'.format(self.model.objective_value))
        
        self.cp = CutPool()
        while continueSeparation:
            iterations = iterations + 1
            valid_cuts = 0
            for d in delta:
                valid_cuts = valid_cuts + self.splitCutsMip(d)
                # print('{} '.format(d), end='')
                
            print(valid_cuts)
            
            if valid_cuts == 0:
                if flag == 1:
                    continueSeparation = False
                else:
                    flag = 1
                    # enrich grid
                    aux = [0]*(2*len(grid)-1)
                    aux_delta = list()
                    for i in range(len(grid)-1):
                        aux[2*i] = grid[i]
                        aux[2*i +1 ] = (grid[i] + grid[i+1])/2
                        aux_delta.append((grid[i] + grid[i+1])/2) 
                    aux[len(aux)-1] = grid[len(grid)-1]
                    grid = aux
                    delta = aux_delta
                    print('new grid: ', grid)
                    print('new delta: ', delta)
            else:
                for c in self.cp.cuts[cuts:]:
                    cuts = cuts + 1
                    # print(c)
                    self.model += c, 'cut_split({})'.format(cuts)
                flag = 0
                
                # self.model.optimize()
                # print('solucao exata: {}'.format(self.model.objective_value))
                self.model.verbose = 0
                self.model.write('{}_split.lp'.format(self.instance.instancename))
                self.model.optimize(relax=True)
                # self.printSolution()
                print('solucao relaxada: {}'.format(self.model.objective_value))
                # input()
            print('cuts = {}, continue = {}, flag = {}'.format(valid_cuts, continueSeparation, flag))
            
            # input()
        print('total', cuts)
        

    def splitCutsMip(self, d):
        print('d', d)
        m = self.m
        m.clear()
        m.verbose = 0
        
        constrs = [c for c in self.model.constrs if c.name.find('cut') == -1]

        u = [m.add_var(var_type=CONTINUOUS, lb=0, name='u({})'.format(constrs[i].name)) for i in range(len(constrs))]
        v = [m.add_var(var_type=CONTINUOUS, lb=0, name='v({})'.format(constrs[i].name)) for i in range(len(constrs))]
        r = [m.add_var(var_type=INTEGER, lb=0, ub=self.model.vars[i].ub, name='r({})'.format(self.model.vars[i].name)) for i in range(len(self.model.vars))]
        r0 = m.add_var(var_type=INTEGER, lb=0, ub=self.instance.K, name='r0') 

        newconstraints = {self.model.vars[i].name: [0]*2 for i in range(len(self.model.vars))}
        newconstraints['b'] = [0]*2
        # input(newconstraints)
        for i in range(len(constrs)):
            mult = 1
            c = constrs[i]
            if c.expr.sense == '<':
                mult = -1
            for (var, coeff) in c.expr.expr.items():
                newconstraints[var.name][0] = newconstraints[var.name][0] + mult * u[i] * coeff
                newconstraints[var.name][1] = newconstraints[var.name][1] + mult * v[i] * coeff

            newconstraints['b'][0] = newconstraints['b'][0] + (-1*c.expr.const) * u[i]
            newconstraints['b'][1] = newconstraints['b'][1] + (-1*c.expr.const) * v[i]

        lin = 0

        for i in range(len(constrs)):
            c = constrs[i]
            if c.expr.sense == '<':
                mult = -1
            soma = 0

            # print('(', end='')
            for (var, coeff) in c.expr.expr.items():
                # print('aj',(coeff*mult), 'x',var.x, '+ ', end='')
                soma = soma + (coeff*mult)*var.x
            
            soma = soma + c.expr.const * mult
            # print('b', c.expr.const,') =',soma,u[i])
            
            lin = lin + soma*u[i]

        for i in range(len(self.model.vars)):
            # print(d,r[i],self.model.vars[i].x)
            lin = lin - d*r[i]*self.model.vars[i].x
            m += newconstraints[self.model.vars[i].name][0] - newconstraints[self.model.vars[i].name][1] - r[i] == 0, 'c35({})'.format(self.model.vars[i].name)
        lin = lin + d*r0
        m.objective = lin

        m += newconstraints['b'][1] - newconstraints['b'][0] + r0 == d-1, 'c37'
        
        # print(lin)
        # input()
        
        # m.write('teste.lp')
        # m.verbose = 1
        start = time()
        m.optimize()
        end = time()
        print('elapsed: {} sec'.format(round(end-start, 4)))
        # for i in range(len(u)):
        #     print(u[i].name, u[i].x)
        #     print(v[i].name, v[i].x)
        # for i in range(len(r)):
        #     print(r[i].name, r[i].x)
        # input('write_aux')
        
        valid_solutions = 0
        # print(m.num_solutions)
        for k in range(m.num_solutions):
            # print(r0.name, r0.xi(k))
            # print('obj', m.objective_values[k])
            if round(m.objective_values[k], 8) < -1e-8:
                
                verificacao = 0
                lin = 0        
                for i in range(len(constrs)):
                    c = constrs[i]
                    if c.expr.sense == '<':
                        mult = -1
                    soma = 0
                    aux = 0
                    # print('(', end='')
                    for (var, coeff) in c.expr.expr.items():
                        # print('aj',(coeff*mult), 'x',var.x, '+ ', end='')
                        soma = soma + (coeff*mult)*var
                        aux = aux + (coeff*mult)*var.x
                    
                    soma = soma + c.expr.const * mult
                    aux = aux + c.expr.const * mult
                    # print('b', c.expr.const,') =',soma,u[i])
                    
                    verificacao = verificacao + aux*u[i].xi(k)
                    lin = lin + soma*u[i].xi(k)

                for i in range(len(self.model.vars)):
                    # print(d,r[i],self.model.vars[i].x)
                    lin = lin - d*r[i].xi(k)*self.model.vars[i]
                    verificacao = verificacao - d*r[i].xi(k)*self.model.vars[i].x
                lin = lin + d*r0.xi(k)
                verificacao = verificacao + d*r0.xi(k)

                # for i in range(len(constrs)):
                #     c = constrs[i]
                #     if c.expr.sense == '<':
                #         mult = -1
                #     soma = 0
                #     for (var, coeff) in c.expr.expr.items():
                #         verificacao = verificacao + (coeff*mult)*var.x
                #         soma = soma + (coeff*mult)*var
                    
                #     soma = soma + c.expr.const
                #     verificacao = verificacao + c.expr.const
                #     lin = lin + soma*u[i].xi(k)
                #     verificacao = verificacao*u[i].xi(k)
                # for i in range(len(self.model.vars)):
                #     lin = lin - d*r[i].xi(k)*self.model.vars[i]
                #     verificacao = verificacao - d*r[i].xi(k)*self.model.vars[i].x
                # lin = lin + d*r0.xi(k)
                # verificacao = verificacao + d*r0.xi(k)
                # print(lin)
                # self.printSolution()
                if self.cp.add(lin >= 0):
                    valid_solutions = valid_solutions + 1
                    # print('obj', m.objective_values[k], 'verificacao', verificacao)
                # self.model.write('{}1_split.lp'.format(self.instance.instancename))
                # self.model.optimize(relax=True)
                # self.printSolution()
                # input('write')

        return valid_solutions



