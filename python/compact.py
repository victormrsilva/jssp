# compact.py
import string

from mip.model import Model, xsum
from mip.constants import INTEGER, BINARY, CONTINUOUS
import sys
from time import process_time, time
import numpy as np
import copy

# from compact_cutpool import Compact_CutPool
from simulatedannealing import SimulatedAnnealing

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
        self.x = 0
        self.y = 0
        self.v = 0
        self.u = 0

    # build the problem with big-M
    def constructProblemM(self):
        self.instance.print()

        self.c = self.model.add_var(var_type=INTEGER, name="C")
        self.x = [[self.model.add_var(var_type=INTEGER, name='x({},{})'.format(j, i)) for i in range(self.instance.m)]
                  for j in range(self.instance.n)]
        self.y = [
            [[self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(j, k, i)) for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]

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

        self.c = self.model.add_var(var_type=INTEGER, lb=0, ub=self.instance.K, name="C")
        self.x = [[self.model.add_var(var_type=INTEGER, lb=self.instance.est[j][i], ub=self.instance.lst[j][i],
                  name='x({},{})'.format(j, i)) for i in range(self.instance.m)] for j in range(self.instance.n)]
        self.y = [
            [[self.model.add_var(var_type=BINARY, name='y({},{},{})'.format(j, k, i)) for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]

        self.v = [
            [[self.model.add_var(var_type=INTEGER, lb=(self.instance.est[k][i] - self.instance.lst[j][i] - self.instance.times[j][i]),
                   ub=(self.instance.lst[k][i] - self.instance.est[j][i] - self.instance.times[j][i]), name='v({},{},{})'.format(j, k, i)) for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]
        self.u = [
            [[self.model.add_var(var_type=INTEGER, lb=-10*self.instance.K, name='u({},{},{})'.format(j, k, i)) for i in range(self.instance.m)]
             for k in range(self.instance.n)] for j in range(self.instance.n)]

        self.model.objective = self.c

        # constraints (2)
        for j in range(self.instance.n):
            for i in range(1, self.instance.m):
                self.model += self.x[j][self.instance.machines[j][i]] - self.x[j][self.instance.machines[j][i - 1]] >= \
                              self.instance.times[j][self.instance.machines[j][i - 1]], 'ord({},{})'.format(j, i)

        # constraints (3-4)
        for j in range(self.instance.n):
            for k in range(j + 1, self.instance.n):
                for i in range(self.instance.m):
                    vjkL = self.instance.est[k][i] - self.instance.lst[j][i] - self.instance.times[j][i]
                    vjkU = self.instance.lst[k][i] - self.instance.est[j][i] - self.instance.times[j][i]
                    vkjL = self.instance.est[j][i] - self.instance.lst[k][i] - self.instance.times[k][i]
                    vkjU = self.instance.lst[j][i] - self.instance.est[k][i] - self.instance.times[k][i]
                    self.model += self.v[j][k][i] - self.x[k][i] + self.x[j][i] == - self.instance.times[j][i], 'MCLinV({},{},{})'.format(j, k, i)
                    self.model += self.v[k][j][i] - self.x[j][i] + self.x[k][i] == - self.instance.times[k][i], 'MCLinV({},{},{})'.format(k, j, i)
                    self.model += self.u[j][k][i] >= 0, 'MCLin1({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] - vjkL * self.y[j][k][i] >= 0, 'MCLin2({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] - vjkU * self.y[j][k][i] - self.v[j][k][i] >= - vjkU, 'MCLin3({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] - vjkU * self.y[j][k][i] <= 0, 'MCLin4({},{},{})'.format(j, k, i)
                    self.model += self.u[j][k][i] - self.v[j][k][i] - vjkL * self.y[j][k][i] <= - vjkL, 'MCLin5({},{},{})'.format(j, k, i)
                    self.model += self.u[k][j][i] >= 0, 'MCLin6({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] + vkjL * self.y[k][j][i] >= vkjL, 'MCLin7({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] + vkjU * self.y[j][k][i] - self.v[k][j][i] >= 0, 'MCLin8({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] - vkjU * self.y[j][k][i] <= vkjU , 'MCLin9({},{},{})'.format(k, j, i)
                    self.model += self.u[k][j][i] + vkjL * self.y[j][k][i] - self.v[k][j][i] <= 0, 'MCLin10({},{},{})'.format(k, j, i)

        # constraints (5)
        for j in range(self.instance.n):
            self.model += self.c - self.x[j][self.instance.machines[j][self.instance.m - 1]] >= self.instance.times[j][
                self.instance.machines[j][self.instance.m - 1]], 'makespan({})'.format(j)
        self.model.write(
            '{}_modelMCLin.lp'.format(self.instance.instancename.translate(str.maketrans('', '', string.punctuation))))


    # cutpool
    def optmizeCuts(self):
        self.model.cuts_generator = Compact_CutPool(self.instance)
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
        # print('Best basic cuts')
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        e_p = [m.add_var(var_type=INTEGER, lb=0, name='e_p({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        y = [m.add_var(var_type=INTEGER, lb=0, name='y({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='e({})'.format(i)) for i in range(self.instance.n)]
        z = m.add_var(var_type=INTEGER, name='E')

        # var = xsum(self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n))

        # var += - xsum(self.instance.times[j][a]*e_p[j] for j in range(self.instance.n))
        # z*xsum(self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n))
        # var += - xsum(self.instance.times[j][a] * self.instance.times[i][a] * xij[i][j] for i in range(self.instance.n) for j in range(i+1,self.instance.n))
        m.objective = xsum(
            self.instance.times[j][a] * self.x[j][a].x * x_aux[j] for j in range(self.instance.n)) - xsum(
            self.instance.times[j][a] * e_p[j] for j in range(self.instance.n)) - xsum(
            self.instance.times[j][a] * self.instance.times[i][a] * xij[j][i] for i in range(self.instance.n) for j in
            range(i + 1, self.instance.n))
        # print(var)

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

        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # print('{}'.format((round(m.objective_value,2))) )
        # for j in range(self.instance.n):
        #     print('x[{}] = {}'.format(j, round(x_aux[j].x,2)))
        #     print('v[{}] = {}'.format(j, round(v[j].x,2)))
        #     print('y[{}] = {}'.format(j, round(y[j].x,2)))
        #     print('o[{}] = {}'.format(j, round(o[j].x,2)))
        # print('z = {}'.format(round(z.x,2)))
        # input()

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
        # print(xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls)
        self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, 'basic_cuts_best' \
                                                                                       '{}({},{})'.format(
            self.iterationsCuts, ''.join(str(i) for i in S), a)
        # input()
        return 1

    def basic_cuts_plus_epsilon_best(self, a, k):
        # print('Best basic cuts plus epsilon. K = {}'.format(k))
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xij = [[m.add_var(var_type=BINARY, lb=0, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)] for j
               in range(self.instance.n)]
        for i in range(self.instance.n):
            for j in range(i + 1, self.instance.n):
                xij[i][j] = xij[j][i]

        # for j in range(self.instance.n):
        #     if j != k:
        #         print('+{}*{}*{}'.format(self.instance.times[j][a], self.x[j][a].x , x_aux[j].name), end='')
        # print()
        # for i in range(self.instance.n):
        #     for j in range(i+1, self.instance.n):
        #         if j != k and i != k:
        #             print('- {}*{}*{}'.format(self.instance.times[j][a], self.instance.times[i][a], xij[i][j].name), end='')
        # print()
        # for j in range(self.instance.n):
        #     if k != j:
        #         print('-{}*{}*{}'.format(self.instance.e[k][a], self.instance.times[j][a], x_aux[j].name), end='')
        # print()
        # for j in range(self.instance.n):
        #     if j != k:
        #         print('- {}*{}*{}*{}'.format(self.y[j][k][a].x, max(self.instance.e[k][a] - self.instance.e[j][a], 0), self.instance.times[j][a], x_aux[j].name), end='')
        # print()
        # for i in range(self.instance.n):
        #     for j in range(self.instance.n):
        #         if j != k and i != k and j != k:
        #             print('- {}*{}*{}*{}'.format(self.y[i][k][a].x, max(self.instance.e[k][a] - self.instance.e[i][a], 0), self.instance.times[j][a], xij[i][j].name), end='')
        # print()

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
        # print(var)
        # input()
        m.objective = var
        # print(var)

        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += x_aux[k] == 0, 'not_in_S'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'

        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # print('{}'.format((round(m.objective_value,2))) )
        # for j in range(self.instance.n):
        #     print('x[{}] = {}'.format(j, round(x_aux[j].x,2)))
        # input()

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

        # print(xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,a)>= ls)
        self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(
            self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,
                                                                                                         a) >= ls, 'basic_cuts_plus_epsilon_best{}({},{},{})'.format(
            self.iterationsCuts, ''.join(str(i) for i in S), a, k)
        # input()
        return 1

    def half_cuts_best(self, a, k):
        # print('Best half cuts. K = {}'.format(k))
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        t = [m.add_var(var_type=INTEGER, lb=0, name='t({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='o({})'.format(i)) for i in range(self.instance.n)]
        e = m.add_var(var_type=INTEGER, name='E')
        C = m.add_var(name='C', lb=-100 * self.instance.K)

        # print(self.x[k][a].x)
        # print(e)
        # for j in range(self.instance.n):
        #     if j != k:
        #         print(' - {}*{}*{}'.format( self.y[j][k][a].x,self.instance.times[j][a],x_aux[j].name,end=''))
        # print()
        var = self.x[k][a].x - e - xsum(
            self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k)

        m.objective = C
        # print(var)

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

        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # print('{}'.format((round(m.objective_value,2))) )
        # for j in range(self.instance.n):
        #     print('x[{}] = {}'.format(j, round(x_aux[j].x,2)))
        #     print('v[{}] = {}'.format(j, round(v[j].x,2)))
        #     print('t[{}] = {}'.format(j, round(t[j].x,2)))
        #     print('o[{}] = {}'.format(j, round(o[j].x,2)))
        # print('e = {}'.format(round(e.x,2)))
        # input()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        self.model += self.x[k][a] - self.E(S, a) - xsum(
            self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) >= 0, 'half_cuts_best{}({},{},{})'.format(
            self.iterationsCuts, k, ''.join(str(i) for i in S), a)
        # input()
        return 1

    def late_job_cuts_best(self, a, k, l):
        # print('Best Late jobs cuts. K = {}. L = {}'.format(k, l))
        m = self.m
        m.clear()
        m.verbose = 0
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        C = m.add_var(name='C', lb=-100 * self.instance.K)

        # print(self.x[k][a].x)
        # print(- self.instance.e[l][a])
        # for j in range(self.instance.n):
        #     if j != k:
        #         print(' - {}*{}*{}'.format( self.y[j][k][a].x,self.instance.times[j][a],x_aux[j].name,end=''))
        # print()
        # for j in range(self.instance.n):
        #     if j != l:
        #         print(' + {}*{}*{}'.format( self.y[j][l][a].x,max(self.instance.e[l][a] - self.instance.e[j][a], 0),x_aux[j].name,end=''))
        # print()
        var = self.x[k][a].x - self.instance.e[l][a] - xsum(
            self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k)
        var += xsum(self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0) * x_aux[j] for j in
                    range(self.instance.n) if j != l)

        m.objective = C
        # print(var)

        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'
        m += C + xsum(self.y[j][k][a].x * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n)
                      if j != k) - xsum(self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0)
                                        * x_aux[j] for j in range(self.instance.n) if j != l) == self.x[k][a].x - \
             self.instance.e[l][a]

        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # print('obj: {}'.format((round(m.objective_value,2))) )
        # for j in range(self.instance.n):
        #     print('x[{}] = {}'.format(j, round(x_aux[j].x,2)))
        # input()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0

        self.model += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
            self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
                      self.instance.e[l][a], 'late_job_cuts_best{}({},{},{},{})'.format(self.iterationsCuts,
                                                                                        ''.join(str(i) for i in S), a,
                                                                                        k, l)
        # input('corte criado')
        return 1

    def two_jobs_cuts_best(self, a):
        # print('Best two jobs cuts')
        m = self.m
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
        # print(var)
        # input()
        m.objective = var

        # print(var)

        for j in range(self.instance.n):
            for i in range(j + 1, self.instance.n):
                m += xij[i][j] <= x_aux[i], 'xij1({},{})'.format(j, i)
                m += xij[i][j] <= x_aux[j], 'xij2({},{})'.format(j, i)
                m += xij[i][j] >= x_aux[j] + x_aux[i] - 1, 'xij3({},{})'.format(j, i)
        m += xsum(x_aux[j] for j in range(self.instance.n)) == 2, 'must_be_2'

        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # print('{}'.format((round(m.objective_value,2))) )
        # for j in range(self.instance.n):
        #     print('x[{}] = {}'.format(j, round(x_aux[j].x,2)))
        # input()

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
        rs = (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a].x + (
                self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a].x
        # left side
        ls = self.instance.times[i][a] * self.instance.times[j][a] + self.instance.e[i][a] * \
             self.instance.times[j][a] + self.instance.e[j][a] * self.instance.times[i][a]

        # apply cut to the model
        self.model += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
            a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
                          a] >= ls, 'two_jobs_cuts{}({},{},{})'.format(self.iterationsCuts, i, j, a)
        return 1

    def clique_cuts(self, S, a):
        s_aux = list(range(len(S)))
        perms = np.asarray(list(permutations(s_aux)))
        perms = [np.asarray(s) for s in perms]
        # print(perms)

        K = [[0] * len(S) for i in range(len(perms))]
        # print(K)

        for i in range(len(K)):
            soma = 0
            for j in range(len(K[i])):
                if j == 0:
                    soma += self.instance.e[S[perms[i][j]]][a]
                else:
                    soma += self.instance.times[S[perms[i][j - 1]]][a]
                for aux in range(j, len(K[i])):
                    K[i][perms[i][aux]] = soma
        # print(K)
        # input('input')
        m = self.m
        m.clear()
        m.verbose = 0
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        m.objective = xsum(i for i in t)

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)

        # m.write('model_clique.lp')
        m.optimize()

        soma = 0
        for i in range(len(S)):
            soma += t[i].x * self.x[S[i]][a].x

        cuts = 0
        # xt >= is a valid inequality. Thus, a violated cut is < 1
        if (1 - soma) > 0.00001:  # < 1:
            cuts = 1
            self.model += xsum(t[i].x * self.x[S[i]][a] for i in range(len(S))) >= 1, 'clique_cuts{}({},{})'.format(
                self.iterationsCuts, ''.join(str(i) for i in S), a)

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
                # print('{}: {}'.format(S, dist))
                dict[s] = dist
            # print(dict)
            # print(sorted(dict.items(), key=lambda l: l[1]))
            i = 0
            # print('Elapsed time: {}'.format(end-start))
            end = time()
            if (end - start) >= timeLimit:
                print('Time limit for enumaration of cliques')
                return 0

            for key, value in sorted(dict.items(), key=lambda l: l[1]):
                end = time()
                if (end - start) >= timeLimit:
                    timeout = True
                    print('Time limit for finding cliques. Cliques found: {}'.format(cliquesFound))
                    return cliquesFound
                S = comb[key]
                # print(S)
                cliquesFound += self.clique_cuts(S, a)
                i += 1
                if i == self.maxCliques:
                    break
            # input(cliquesFound)
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
        xij = [[m.add_var(var_type=BINARY, name='xij({},{})'.format(i, j)) for i in range(self.instance.n)]
               for j in range(self.instance.n)]

        var = - xsum(xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if i != j)
        # print(var)
        # input()
        m.objective = var
        # print(var)

        for i in range(self.instance.n):
            m += xsum(xij[i][j] for j in range(self.instance.n) if j != i) - xsum(xij[j][i]
                                                                                  for j in range(self.instance.n) if
                                                                                  j != i) == 0, 'flow({})'.format(i)
        m += xsum(
            self.y[i][j][a].x * xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i
            ) - xsum(
            xij[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if j != i) >= -1, 'triangle_cut'
        # m.write('best_model.lp')
        # input('modelo pronto')
        m.optimize()
        # if m.objective_value < 0:
        #     m.write('best_model.lp')
        #     print('{}'.format((round(m.objective_value,2))) )
        #     for i in range(self.instance.n):
        #         for j in range(self.instance.n):
        #             print('x[{}][{}] = {} - y[{},{}] = {}'.format(i, j, round(xij[i][j].x, 2), i, j, round(self.y[i][j][a].x, 2)))
        #     input()

        if m.objective_value > -0.000000001:
            return 0

        S = []
        soma = 0
        for i in range(self.instance.n):
            for j in range(self.instance.n):
                if xij[i][j].x > 0.99999999:
                    S.append(self.y[i][j][a])
                    soma += self.y[i][j][a].x * xij[i][j].x

        if len(S) <= 1 or soma - (len(S) - 1) <= 0.000000001:
            # print('Soma: {}. Len: {}'.format(soma, len(S)))
            # input()
            return 0

        self.model += xsum(i for i in S) <= len(S) - 1, 'triangle_cuts_best{}({})'.format(self.iterationsCuts, a)
        # for i in S:
        #     print('{} = {}; '.format(i.name, i.x), end='')
        #
        # print('Soma: {}. Len: {}'.format(soma, len(S)))

        # input()
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
                print('x({},{}) = {} '.format(j, i, self.x[j][i].x), end='')
            print()
