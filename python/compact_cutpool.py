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

            result = self.basic_cuts_best(a)
            if isinstance(result,int) == False:
                print(result)
                cp.add(result)

            for k in range(self.instance.n):
                result = self.basic_cuts_plus_epsilon_best(a,k)
                if isinstance(result,int) == False:
                    print(result)
                    cp.add(result)

                result = self.half_cuts_best(a,k)
                if isinstance(result,int) == False:
                    print(result)
                    cp.add(result)

                for l in range(self.instance.n):
                    result = self.late_job_cuts_best(a,k,l)
                    if isinstance(result,int) == False:
                        print(result)
                        cp.add(result)

            # for sizeS in range(1, self.instance.n):
            #     # generate all combinations of jobs with size sizeS = {2,...,n}
            #     # used for cuts where the order of jobs doesn't matter
            #     comb = combinations(list(range(0, self.instance.n)), sizeS + 1)
            #     for S in list(comb):
            #         result = self.basic_cuts(S,a)
            #         if isinstance(result,int) == False:
            #             print(result)
            #             cp.add(result)

            #         result = self.basic_cuts_reverse(S,a)
            #         if isinstance(result,int) == False:
            #             print(result)
            #             cp.add(result)

            #         for k in range(self.instance.n):
            #             result = self.basic_cuts_plus_epsilon(S,a,k)
            #             if isinstance(result,int) == False:
            #                 print(result)
            #                 cp.add(result)

            #         for k in S:
            #             for l in range(self.instance.n):
            #                 result = self.late_job_cuts(S,a,k,l)
            #                 if isinstance(result,int) == False:
            #                     print(result)
            #                     cp.add(result)

            #             result = self.half_cuts(S,a,k)
            #             if isinstance(result,int) == False:
            #                 print(result)
            #                 cp.add(result)

            #         if len(S) <= 5:
            #             result = self.clique_cuts(S, a)
            #             if isinstance(result,int) == False:
            #                 print(result)
            #                 cp.add(result)
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

    def basic_cuts_best(self, a):
        m = Model('basic_cut_best')
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        e_p = [m.add_var(var_type=INTEGER, lb=0, name='e_p({})'.format(i)) for i in range(self.instance.n)]
        xi_xj = [[m.add_var(var_type=BINARY, lb=0, name='xi_xj({},{})'.format(i,j)) for i in range(self.instance.n)] for j in range(self.instance.n)]
        y = [m.add_var(var_type=INTEGER, lb=0, name='y({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='e({})'.format(i)) for i in range(self.instance.n)]
        z = m.add_var(var_type=INTEGER, name='E')

        m.objective = xsum(self.instance.times[j][a] * self.x[j][a][1] * x_aux[j] for j in range(self.instance.n)) - xsum(self.instance.times[j][a]*e_p[j] for j in range(self.instance.n)) - xsum(self.instance.times[j][a] * self.instance.times[i][a] * xi_xj[j][i] for i in range(self.instance.n) for j in range(i+1,self.instance.n))

        for j in range(self.instance.n):
            m += v[j] - self.instance.e[j][a]*x_aux[j] + self.instance.K * x_aux[j] == self.instance.K, 'eq26({})'.format(j)
            m += z - v[j] <= 0, 'eq27({})'.format(j)
            m += z - y[j]>= 0, 'eq28({})'.format(j)
            m += y[j] - self.instance.K * o[j] <= 0, 'eq29({})'.format(j)
            m += y[j] - v[j] <= 0, 'eq30({})'.format(j)
            m += y[j] - v[j] - self.instance.K*o[j] >= - self.instance.K, 'eq31({})'.format(j)
            #z*x(j)
            m += e_p[j] - self.instance.K * x_aux[j]<= 0, 'e_s1({})'.format(j)
            m += e_p[j] - z <= 0, 'e_s2({})'.format(j)
            m += e_p[j] - z - self.instance.K * x_aux[j] >= -self.instance.K , 'e_s3({})'.format(j)
            # x[i]*x[j]
            for i in range(j+1, self.instance.n):
                m += xi_xj[i][j] <= x_aux[i], 'xi_xj1({},{})'.format(j,i)
                m += xi_xj[i][j] <= x_aux[j], 'xi_xj2({},{})'.format(j,i)
                m += xi_xj[i][j] >= x_aux[j]+x_aux[i] - 1, 'xi_xj3({},{})'.format(j,i)
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
        
        return xsum(self.instance.times[j][a] * self.x[j][a][0] for j in S) >= ls

    def basic_cuts_plus_epsilon_best(self, a, k):
        m = Model('basic_cut_plus_epsilon_best')
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        xi_xj = [[m.add_var(var_type=BINARY, lb=0, name='xi_xj({},{})'.format(i,j)) for i in range(self.instance.n)] for j in range(self.instance.n)]
        for i in range(self.instance.n):
            for j in range(i+1, self.instance.n):
                xi_xj[i][j] = xi_xj[j][i]

        var = xsum(self.instance.times[j][a] * self.x[j][a][1] * x_aux[j] for j in range(self.instance.n)) # part 1 of cut
        var += - xsum(self.instance.times[j][a] * self.instance.times[i][a] * xi_xj[i][j] for i in range(self.instance.n) for j in range(i+1,self.instance.n)) # part 2 of cut
        var += - self.instance.e[k][a]*xsum(self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k) # part 3 of cut
        var += xsum(self.y[j][k][a][1] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k) # part 4 of cut if i == j (just expand the multiplication)
        var += xsum(self.y[i][k][a][1] * max(self.instance.e[k][a] - self.instance.e[i][a], 0) * self.instance.times[j][a]* xi_xj[i][j] for i in range(self.instance.n) for j in range(self.instance.n) if i != k and j != k and i != j) # part 4 of cut if i != j (just expand the multiplication)
        m.objective = var

        for j in range(self.instance.n):
            for i in range(j+1, self.instance.n):
                m += xi_xj[i][j] <= x_aux[i], 'xi_xj1({},{})'.format(j,i)
                m += xi_xj[i][j] <= x_aux[j], 'xi_xj2({},{})'.format(j,i)
                m += xi_xj[i][j] >= x_aux[j]+x_aux[i] - 1, 'xi_xj3({},{})'.format(j,i)
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
        
        return xsum(self.instance.times[j][a] * self.x[j][a][0] for j in S) + xsum(self.y[j][k][a][0] * max(self.instance.e[k][a] - self.instance.e[j][a], 0) for j in S) * self.p(S,a) >= ls

    def half_cuts_best(self,a,k):
        m = Model('basic_cut_best')
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        v = [m.add_var(var_type=INTEGER, lb=0, name='v({})'.format(i)) for i in range(self.instance.n)]
        t = [m.add_var(var_type=INTEGER, lb=0, name='t({})'.format(i)) for i in range(self.instance.n)]
        o = [m.add_var(var_type=BINARY, lb=0, name='o({})'.format(i)) for i in range(self.instance.n)]
        e = m.add_var(var_type=INTEGER, name='E')
        C = m.add_var(name='C', lb=-100*self.instance.K)
        
        var = self.x[k][a][1] - e - xsum(self.y[j][k][a][1]*self.instance.times[j][a]*x_aux[j] for j in range(self.instance.n) if j != k)

        m.objective = C

        for j in range(self.instance.n):
            m += v[j] - self.instance.e[j][a]*x_aux[j] + self.instance.K * x_aux[j] == self.instance.K, 'eq26({})'.format(j)
            m += e - v[j] <= 0, 'eq27({})'.format(j)
            m += e - t[j]>= 0, 'eq28({})'.format(j)
            m += t[j] - self.instance.K * o[j] <= 0, 'eq29({})'.format(j)
            m += t[j] - v[j] <= 0, 'eq30({})'.format(j)
            m += t[j] - v[j] - self.instance.K*o[j] >= - self.instance.K, 'eq31({})'.format(j)
        m += xsum(o[j] for j in range(self.instance.n)) == 1, 'eq32'
        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'
        m += C + e + xsum(self.y[j][k][a][1]*self.instance.times[j][a]*x_aux[j] for j in range(self.instance.n) if j != k) == self.x[k][a][1]

        m.optimize()

        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0
        
        return self.x[k][a][0] - self.E(S,a) - xsum(self.y[j][k][a][0]*self.instance.times[j][a] for j in S if j != k) >= 0

    def late_job_cuts_best(self,a,k,l):
        m = Model('basic_cut_best')
        x_aux = [m.add_var(var_type=BINARY, lb=0, name='x({})'.format(i)) for i in range(self.instance.n)]
        C = m.add_var(name='C', lb=-100*self.instance.K)
        
        var = self.x[k][a][1] - self.instance.e[l][a] - xsum( self.y[j][k][a][1] * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k)
        var += xsum(self.y[j][l][a][1] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) * x_aux[j] for j in range(self.instance.n) if j != l)

        m.objective = C

        m += xsum(x_aux[j] for j in range(self.instance.n)) >= 2, 'minimum_jobs'
        m += x_aux[k] == 1, 'k_in_S'
        m += C + xsum( self.y[j][k][a][1] * self.instance.times[j][a] * x_aux[j] for j in range(self.instance.n) if j != k) -  xsum(self.y[j][l][a][1] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) * x_aux[j] for j in range(self.instance.n) if j != l) == self.x[k][a][1] - self.instance.e[l][a]

        m.optimize()
        if m.objective_value > -0.0001:
            return 0

        S = []
        for j in range(self.instance.n):
            if x_aux[j].x > 0.99999:
                S.append(j)

        if len(S) <= 1:
            return 0
        
        return self.x[k][a][0] - xsum(self.y[j][k][a][0] * self.instance.times[j][a] for j in S if j != k) + xsum(self.y[j][l][a][0] * max(self.instance.e[l][a] - self.instance.e[j][a], 0) for j in S if j != l) >= self.instance.e[l][a]


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
