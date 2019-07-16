# compact_cutpool.py
from mip.model import Model, xsum
from mip.constants import INTEGER, BINARY, CONTINUOUS
from mip.callbacks import CutsGenerator, CutPool
import sys
from time import process_time
import numpy as np

from itertools import permutations, combinations


class Compact_CutPool(CutsGenerator):
    def __init__(self, instance):
        self.instance = instance

    def generate_cuts(self, model: Model):
        self.x = [[(0,0) for i in range(self.instance.m)] for j in range(self.instance.n)]
        self.y = [[[(0,0) for i in range(self.instance.m)] for k in range(self.instance.n)] for j in range(self.instance.n)]

        # recreate the set of vars used in the model. if position 1 have value = 0 and not a var object
        # the variable isn't used in the model
        for v in model.vars:
            if v.name.startswith('C'):
                self.c = (v,v.x)
            if v.name.startswith('x('):
                j = int(v.name.split('(')[1].split(',')[0]) - 1
                i = int(v.name.split(')')[0].split(',')[1]) - 1
                self.x[j][i] = (v, v.x)
            if v.name.startswith('y('):
                j = int(v.name.split('(')[1].split(',')[0]) - 1
                k = int(v.name.split('(')[1].split(',')[1]) - 1
                i = int(v.name.split(')')[0].split(',')[1]) - 1
                self.y[j][k][i] = (v, v.x)
        # print('{} = {}'.format(self.c[0].name, self.c[1]))
        # for j in range(self.instance.n):
        #     for i in range(self.instance.m):
        #         if isinstance(self.x[j][i][0], int) == False:
        #             print('{} = {}'.format(self.x[j][i][0].name, self.x[j][i][1]))
        # for j in range(self.instance.n):
        #     for k in range(self.instance.n):
        #         for i in range(self.instance.m):
        #             if isinstance(self.y[j][k][i][0], int) == False:
        #                 print('{} = {}'.format(self.y[j][k][i][0].name, self.y[j][k][i][1]))
        #
        # input()
        cp = CutPool()
        for a in range(self.instance.m):
            # all possible combinations with 2 jobs. order doesnt' matter
            comb = combinations(list(range(0,self.instance.n)),2)
            for jobs in list(comb):
                result = self.two_jobs_cuts(jobs[0],jobs[1],a)
                if isinstance(result,int) == False:
                    print(result)
                    cp.add(result)
            # all possible combinations with 3 josb. the order matters
            perm = permutations(list(range(0, self.instance.n)), 3)
            for jobs in list(perm):
                result = self.triangle_cuts(jobs[0],jobs[1], jobs[2], a)
                if isinstance(result,int) == False:
                    print(result)
                    cp.add(result)

            for sizeS in range(1, self.instance.n):
                # generate all combinations of jobs with size sizeS = {2,...,n}
                # used for cuts where the order of jobs doesn't matter
                comb = combinations(list(range(0, self.instance.n)), sizeS + 1)
                for S in list(comb):
                    result = self.basic_cuts(S,a)
                    if isinstance(result,int) == False:
                        print(result)
                        cp.add(result)

                    result = self.basic_cuts_reverse(S,a)
                    if isinstance(result,int) == False:
                        print(result)
                        cp.add(result)

                    for k in range(self.instance.n):
                        result = self.basic_cuts_plus_epsilon(S,a,k)
                        if isinstance(result,int) == False:
                            print(result)
                            cp.add(result)

                    for k in S:
                        for l in range(self.instance.n):
                            result = self.late_job_cuts(S,a,k,l)
                            if isinstance(result,int) == False:
                                print(result)
                                cp.add(result)

                        result = self.half_cuts(S,a,k)
                        if isinstance(result,int) == False:
                            print(result)
                            cp.add(result)

                    if len(S) <= 5:
                        result = self.clique_cuts(S, a)
                        if isinstance(result,int) == False:
                            print(result)
                            cp.add(result)
        for cut in cp.cuts:
            model.add_cut(cut)
        return

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
            rs += self.instance.times[j][a] * self.x[j][a][1]

        # left side
        ls = self.E(S, a) * self.p(S, a)
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]
        # check violated cut
        if (ls - rs) > 0.00001: # if rs < ls:
            return xsum(self.instance.times[j][a] * self.x[j][a][0] for j in S) >= ls

        return 0

    def basic_cuts_reverse(self, S, a):
        # reverse
        # right side
        rs = 0
        for j in S:
            rs += self.instance.times[j][a] * (self.c[1] - self.x[j][a][1])
        # left side
        ls = self.F(S, a) * self.p(S, a)
        for j in S:
            ls += self.instance.times[j][a] * self.instance.times[j][a]
        for i in range(len(S)):
            for j in range(i + 1, len(S)):
                ls += self.instance.times[S[i]][a] * self.instance.times[S[j]][a]

        # check violated cut
        if (ls - rs) > 0.00001: # if rs < ls:
            return xsum(self.instance.times[j][a] * (self.c[0] - self.x[j][a][0]) for j in S) >= ls

        return 0

    # two-job cuts
    def two_jobs_cuts(self, i, j, a):
        # right side
        rs = (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a][1] + (
                    self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a][1]
        # left side
        ls = self.instance.times[i][a] * self.instance.times[j][a] + self.instance.e[i][a] * \
             self.instance.times[j][a] + self.instance.e[j][a] * self.instance.times[i][a]

        # check violated cut
        if (ls - rs) > 0.00001: # if rs < ls:
            return (self.instance.times[i][a] + self.instance.e[i][a] - self.instance.e[j][a]) * self.x[i][a][0] + (self.instance.times[j][a] + self.instance.e[j][a] - self.instance.e[i][a]) * self.x[j][a][0] >= ls

        return 0

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

        m = Model('clique_cuts')
        t = [m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S]
        m.objective = xsum(i for i in t)

        for i in range(len(K)):
            m += xsum(K[i][j] * t[j] for j in range(len(K[i]))) >= 1, 'K({})'.format(i)

        m.optimize()

        soma = 0
        for i in range(len(S)):
            soma += t[i].x * self.x[S[i]][a][1]

        # xt >= is a valid inequality. Thus, a violated cut is < 1
        if (1 - soma) > 0.00001: # < 1:
            return xsum(t[i].x * self.x[S[i]][a][0] for i in range(len(S))) >= 1

        return 0

    # triangle cuts
    def triangle_cuts(self, i, j, k, a):
        if isinstance(self.y[i][j][a][0],int) or isinstance(self.y[j][k][a][0],int) or isinstance(self.y[k][i][a][0],int):
            return 0
        soma = self.y[i][j][a][1] + self.y[j][k][a][1] + self.y[k][i][a][1]
        if soma > 2:
            return self.y[i][j][a][0] + self.y[j][k][a][0] + self.y[k][i][a][0] <= 2
        return 0

    def basic_cuts_plus_epsilon(self, S, a, k):
        # right side
        rs = 0
        for j in S:
            rs += self.instance.times[j][a] * self.x[j][a][1]
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
            if (isinstance(self.y[j][k][a][0],int) == False):
                soma2 += self.y[j][k][a][1] * max(self.instance.e[k][a] - self.instance.e[j][a], 0)
            else:
                return 0
        ls -= soma2 * self.p(S, a)
        cuts = 0

        if (ls - rs) > 0.00001: # if rs < ls:
            return xsum(self.instance.times[j][a] * self.x[j][a][0] for j in S) + xsum(self.y[j][k][a][0] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,a) >= soma1

        return 0

    def half_cuts(self, S, a, k):
        # left side
        ls = self.E(S, a)
        for j in S:
            if isinstance(self.y[j][k][a][0], int) == True:
                return 0
            if j != k:
                ls += self.y[j][k][a][1] * self.instance.times[j][a]
        if (ls - self.x[k][a][1]) > 0.00001: #self.x[k][a][1] < ls:
            return self.x[k][a][0] - xsum(self.y[j][k][a][0]* self.instance.times[j][a] for j in S if j != k) >= self.E(S,a)

        return 0

    def late_job_cuts(self, S, a, k, l):
        # left side
        ls = self.instance.e[l][a]
        for j in S:
            if isinstance(self.y[j][k][a][0], int) == True:
                return 0
            if j != k:
                ls += self.y[j][k][a][1] * self.instance.times[j][a]
        for j in S:
            if isinstance(self.y[j][l][a][0], int) == True:
                return 0
            if j != l:
                ls -= self.y[j][l][a][1] * max(self.instance.e[l][a] - self.instance.e[j][a], 0)
        if (ls - self.x[k][a][1]) > 0.00001: # self.x[k][a][1] < ls:
            return self.x[k][a][0] - xsum(self.y[j][k][a][0] * self.instance.times[j][a] for j in S if j != k ) + xsum(
                self.y[j][l][a][0] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l ) >= self.instance.e[l][a]

        return 0
