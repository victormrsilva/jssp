# compact.py
from mip.model import Model, xsum
from mip.constants import INTEGER, BINARY, CONTINUOUS
import sys
from time import process_time
import numpy as np

from itertools import permutations, combinations


class Compact:
    def __init__(self, instance):
        self.instance = instance

    # build the problem
    def constructProblem(self):
        self.instance.print()
        # S = [0, 2]
        # print("E({0,2},1) = ", self.E(S, 1))
        # print("F({0,2},1) = ", self.F(S, 1))
        # print("p({0,2},1) = ", self.p(S, 1))
        # S = [0, 1, 2]
        # print("E({0,1,2},0) = ", self.E(S, 0))
        # print("F({0,1,2},0) = ", self.F(S, 0))
        # print("p({0,1,2},0) = ", self.p(S, 0))

        # for aux in range(1,self.instance.m):
        #     comb = combinations(list(range(0,self.instance.n)), aux+1)
        #     for S in list(comb):
        #         print(S)

        self.model = Model('compact')

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
        newConstraints = True
        self.model.write('model.lp')

        self.model.relax()
        self.iterationsCuts = 0
        while newConstraints:
            self.iterationsCuts += 1
            newConstraints = False
            self.model.optimize()
            print('objective value : {}'.format(self.model.objective_value))
            for j in range(self.instance.n):
                for i in range(self.instance.m):
                    print('x({},{}) = {}\t'.format(j, i, round(self.x[j][i].x, 2)), end='')
                print()
            self.model.write('teste.lp')
            hasCuts = 0
            for a in range(self.instance.m):
                triangle_cuts = 0
                basic_cuts_epsilon = 0
                basic_cuts = 0
                late_jobs_cuts = 0
                half_cuts = 0
                two_job_cuts = 0
                clique_cuts = 0
                # all possible combinations with 2 jobs. order doesnt' matter
                comb = combinations(list(range(0,self.instance.n)),2)
                for jobs in list(comb):
                    two_job_cuts += self.two_jobs_cuts(jobs[0],jobs[1],a)
                # all possible combinations with 3 josb. the order matters
                perm = permutations(list(range(0, self.instance.n)), 3)
                for jobs in list(perm):
                    triangle_cuts += self.triangle_cuts(jobs[0],jobs[1], jobs[2], a)

                for sizeS in range(1, self.instance.n):
                    # generate all combinations of jobs with size sizeS = {2,...,n}
                    # used for cuts where the order of jobs doesn't matter
                    comb = combinations(list(range(0, self.instance.n)), sizeS + 1)
                    for S in list(comb):
                        basic_cuts += self.basic_cuts(S,a)

                        for k in range(self.instance.n):
                            # if k not in S:
                            basic_cuts_epsilon += self.basic_cuts_plus_epsilon(S,a,k)

                        for k in S:
                            for l in range(self.instance.n):
                                late_jobs_cuts += self.late_job_cuts(S,a,k,l)
                            half_cuts += self.half_cuts(S,a,k)
                minS = min(self.instance.n,5)
                for sizeS in range(1, minS):
                    # generate all combinations of jobs with size sizeS = {2,...,min(n,5)}
                    # used for cuts where the order of jobs doesn't matter
                    perm = permutations(list(range(0, self.instance.n)), sizeS + 1)
                    for S in list(perm):
                        clique_cuts += self.clique_cuts(S, a)
                print('For machine {} were found {} basic cuts, {} two jobs cuts, {} clique cuts, {} triangle cuts, {} basic cuts epsilon, {} half cuts and {} late jobs cuts'.format(a, basic_cuts, two_job_cuts, clique_cuts, triangle_cuts, basic_cuts_epsilon, half_cuts, late_jobs_cuts))
                hasCuts += basic_cuts + clique_cuts + basic_cuts_epsilon + half_cuts + late_jobs_cuts + two_job_cuts + triangle_cuts
                self.model.write('teste.lp')
                # input()

            if hasCuts > 0:
                newConstraints = True
            self.model.relax()
            self.optimize()

    # optimize model
    def optimize(self):
        self.model.optimize()

        # printing results
        print("C: ", self.c.x)
        for j in range(self.instance.n):
            for i in range(self.instance.m):
                print('x({},{}) = {} '.format(j + 1, i + 1, self.x[j][i].x), end='')
            print()

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

    def basic_cuts(self, S, a):
        cuts = 0

        # right side
        rs = 0
        for j in S:
            # print('+ {}*{} '.format(self.instance.times[j][a], round(self.x[j][a].x, 2)), end='')
            rs += self.instance.times[j][a] * self.x[j][a].x

        # print(' >= ', end='')
        # left side
        ls = self.E(S, a) * self.p(S, a)
        # print(' {}*{} '.format(self.E(S, a), self.p(S, a)), end='')
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                # print('+ {}*{} '.format(self.instance.times[S[i]][a], self.instance.times[S[j]][a]), end='')
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        # print()
        # check violated cut
        if (ls - rs) > 0.00001: # if rs < ls:
            print('{} >= {}'.format(rs, ls))
            cuts += 1
            # print(S)
            self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) >= ls, 'basic_cuts{}({},{})'.format(self.iterationsCuts,
                ''.join(str(i) for i in S), a)

        # reverse
        # right side
        rs = 0
        for j in S:
            # print('+ {}*({} - {}) '.format(self.instance.times[j][a], round(self.c.x, 2), round(self.x[j][a].x, 2)),end='')
            rs += self.instance.times[j][a] * (self.c.x - self.x[j][a].x)
        # print(' >= ', end='')
        # left side
        ls = self.F(S, a) * self.p(S, a)
        # print(' {}*{} '.format(self.F(S, a), self.p(S, a)), end='')
        for j in S:
            # print('+ {}^2 '.format(self.instance.times[j][a]), end='')
            ls += self.instance.times[j][a] * self.instance.times[j][a]
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                # print('+ {}*{} '.format(self.instance.times[S[i]][a], self.instance.times[S[j]][a]), end='')
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        # print()
        # print('{} >= {}'.format(rs, ls))
        # check violated cut

        if (ls - rs) > 0.00001: # if rs < ls:
            print('reverse: {} >= {}'.format(rs, ls))
            cuts += 1
            # print(S)
            self.model += xsum(self.instance.times[j][a] * (self.c - self.x[j][a]) for j in
                               S) >= ls, 'basic_cuts_reverse{}({},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a)

        return cuts

    # two-job cuts
    def two_jobs_cuts(self, i, j, a):
        cuts = 0

        # right side
        rs = (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a].x + (
                    self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a].x
        # left side
        ls = self.instance.times[i][a] * self.instance.times[j][a] + self.instance.e[i][a] * \
             self.instance.times[j][a] + self.instance.e[j][a] * self.instance.times[i][a]

        # print('{} >= {}'.format(rs, ls))
        # input()
        # check violated cut
        if (ls - rs) > 0.00001: # if rs < ls:
            cuts = 1
            self.model += (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][
                a] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][
                              a] >= ls, 'two_jobs_cuts{}({},{},{})'.format(self.iterationsCuts,i, j, a)

        return cuts

    def clique_cuts(self, S, a):
        s_aux = list(range(len(S)))
        perms = np.asarray(list(permutations(s_aux)))
        perms = [np.asarray(s) for s in perms]

        K = [[0] * len(S) for i in range(len(perms))]
        # print(K)
        # print(perms)

        for i in range(len(K)):
            soma = 0
            for j in range(len(K[i])):
                if j == 0:
                    soma += self.instance.e[S[perms[i][j]]][a]
                else:
                    soma += self.instance.times[S[perms[i][j - 1]]][a]
                for aux in range(j, len(K[i])):
                    K[i][perms[i][aux]] = soma

        #print(K)

        m = Model('clique_cuts')
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        m.objective = xsum(i for i in t)

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)

        m.write('model_clique.lp')
        m.optimize()

        soma = 0
        for i in range(len(S)):
            soma += t[i].x * self.x[S[i]][a].x

        cuts = 0
        # xt >= is a valid inequality. Thus, a violated cut is < 1
        if (1 - soma) > 0.00001: # < 1:
            cuts = 1
            # for i in range(len(S)):
            #     print('a({}) = {}\t'.format(t[i].name, t[i].x), end='')
            self.model += xsum(t[i].x * self.x[S[i]][a] for i in range(len(S))) >= 1, 'clique_cuts{}({},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a)
            # input('encontrado')
        return cuts

    # triangle cuts
    def triangle_cuts(self, i, j, k, a):
        soma = self.y[i][j][a].x + self.y[j][k][a].x + self.y[k][i][a].x
        cuts = 0
        if soma > 2:
            cuts = 1
            self.model += self.y[i][j][a] + self.y[j][k][a] + self.y[k][i][a] <= 2, 'triangle_cuts{}({},{},{},{})'.format(self.iterationsCuts,i, j, k, a)
        return cuts

    def basic_cuts_plus_epsilon(self, S, a, k):
        # right side
        rs = 0
        for j in S:
            # print('+ {}*{} '.format(self.instance.times[j][a], round(self.x[j][a].x, 2)), end='')
            rs += self.instance.times[j][a] * self.x[j][a].x
        # print(' >= ', end='')
        # left side
        ls = 0
        soma1 = 0
        for j in range(len(S)):
            for i in range(j + 1, len(S)):
                # print('+ {}*{} '.format(self.instance.times[i][a], self.instance.times[j][a]), end='')
                soma1 += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        # print(' {}*{} '.format(self.instance.e[k][a], self.p(S, a)), end='')
        soma1 += self.instance.e[k][a] * self.p(S, a)

        ls += soma1

        soma2 = 0
        # print(' - {}('.format(self.p(S,a)), end='')
        for j in S:
            # print('+ {} * max({} - {}, 0)'.format(self.y[j][k][a].x, self.instance.e[k][a], self.instance.e[j][a]), end='')
            soma2 += self.y[j][k][a].x * max(self.instance.e[k][a] - self.instance.e[j][a], 0)
        ls -= soma2 * self.p(S, a)
        # print(')')
        # print ('{} >= {}'.format(rs, ls))
        cuts = 0

        if (ls - rs) > 0.00001: # if rs < ls:
            cuts = 1
            self.model += xsum(self.instance.times[j][a] * self.x[j][a] for j in S) + xsum(self.y[j][k][a] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,a) >= soma1, 'basic_cuts_epsilon{}({},{},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a, k)

        return cuts

    def half_cuts(self, S, a, k):
        # print('{} >= '.format(self.x[k][a].x), end='')
        # left side
        ls = self.E(S, a)
        # print('{} '.format(self.E(S,a)), end='')
        for j in S:
            if j != k:
                # print(' + {}*{}'.format(self.y[j][k][a].x,self. instance.times[j][a]), end='')
                ls += self.y[j][k][a].x * self.instance.times[j][a]
        # print()
        # print('{} >= {}'.format(self.x[k][a].x, ls))
        cuts = 0
        if (ls - self.x[k][a].x) > 0.00001: #self.x[k][a].x < ls:
            cuts = 1
            self.model += self.x[k][a] - xsum(self.y[j][k][a]* self.instance.times[j][a] for j in S if j != k) >= self.E(S,a), 'half_cuts{}({},{},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), k, a)
        return cuts

    def late_job_cuts(self, S, a, k, l):
        # left side
        # print('{} >= '.format(self.x[k][a].x), end='')
        # ls = 0
        ls = self.instance.e[l][a]
        # print('{}'.format(self.instance.e[l][a]), end='')
        for j in S:
            if j != k:
                # print(' + {} * {}'.format(self.y[j][k][a].x, self.instance.times[j][a]), end='')
                ls += self.y[j][k][a].x * self.instance.times[j][a]
        for j in S:
            if j != l:
                # print(' - {} * max({} - {},0)'.format(self.y[j][l][a].x, self.instance.e[l][a], self.instance.e[j][a]), end='')
                ls -= self.y[j][l][a].x * max(self.instance.e[l][a] - self.instance.e[j][a], 0)
        # print()
        cuts = 0
        # print('{} >= {}'.format(round(self.x[k][a].x,2), ls))
        if (ls - self.x[k][a].x) > 0.00001: # self.x[k][a].x < ls:
            cuts = 1
            self.model += self.x[k][a] - xsum(self.y[j][k][a] * self.instance.times[j][a] for j in S if j != k) + xsum(
                self.y[j][l][a] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= \
                          self.instance.e[l][a], 'late_job_cuts{}({},{},{},{})'.format(self.iterationsCuts,''.join(str(i) for i in S), a, k,
                                                                                     l)
        return cuts
