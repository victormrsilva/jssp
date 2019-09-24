import numpy as np
import random
import math


class SimulatedAnnealing():

    def __init__(self, temp_i, temp_f, alpha, maxIter, maxSuc, instance, a, x, y):
        """
            temp_i: initial temperature
            temp_f: final temperature
            alpha: alpha
            maxIter: maximum iterations
            maxSuc: maximum success without improvement
            maxPer: maximum Perturbation
            instance: instance of jobshop
            a: machine which the cuts are verified
            x: values of x of the solution
            y: values of y of the solution
        """
        self.temp_i = temp_i
        self.temp_f = temp_f
        self.alpha = alpha
        self.maxIter = maxIter
        self.maxSuc = maxSuc
        self.instance = instance
        self.a = a
        self.x = x
        self.y = y

    def annealing(self, cut):
        """
            si: initial solution
            cut:
                1 for basic cut
                2 for basic cut with epsilon
                3 for ...
                TODO
        """
        s = self.newSolution(cut)  # random solution
        f = self.fo(s, cut)
        t = self.temp_i
        i = 0
        success = 0
        print('Iteration {}, f = {}, t = {}, s = {}'.format(i, f, t, s))
        while i < self.maxIter or success < self.maxSuc:
            s_new, change = self.neighbor(s, cut, f)
            f_new = self.fo(s_new, cut)

            delta = f_new - f

            if f_new < f or np.exp(- (f_new - f) / t) >= random.random():
                s[change] = (s[change] + 1) % 2
                f = f_new
                success = 0
            success = success + 1
            i = i + 1
            t = t * self.alpha
            print('Iteration {}, success = {}, f = {}, f_new = {}, t = {}, s = {}, s_new = {}, change = {}'.format(i,
                        success, f, f_new, t, s, s_new, change))
            input()

    def neighbor(self, s, cut, f):
        change = random.randint(0, self.instance.n - 1)
        aux = s.copy()
        aux[change] = (aux[change] + 1) % 2
        return aux, change

    def fo(self, s, cut):
        fo = 0
        if cut == 1:
            rs = sum(self.instance.times[i][self.a] * self.x[i][self.a].x * s[i] for i in range(len(s)))
            ls = min(self.instance.times[i][self.a] * s[i] for i in range(len(s)) if s[i] != 0) * sum(
                self.instance.times[i][self.a] * s[i] for i in range(len(s))) + sum(
                self.instance.times[i][self.a] * s[i] * self.instance.times[j][self.a] * s[j] for i in range(
                    len(s)) for j in range(i, len(s)))
            fo = rs - ls

        return fo

    def new_fo(self, s, j_change, cut, f):
        return 0

    def newSolution(self, cut):
        if cut == 1:
            return [random.randint(0, 1) for _ in range(self.instance.n)]
