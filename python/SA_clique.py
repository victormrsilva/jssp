import numpy as np
from mip.model import Model, xsum
from mip.entities import Var
from mip.constants import INTEGER, BINARY, CONTINUOUS, OptimizationStatus
from itertools import permutations, combinations
import time

class Clique:
    def __init__(self, x, p, est, maxsteps):
        self.qtd = len(x)
        self.p = p
        self.est = est
        self.x = x
        self.maxsteps = maxsteps
        self.debug = False
        self.max_solutions = 50
        self.m = Model()
        self.accepted = {}
        self.max_exact = 512
        self.exact = 0
        self.minimum = 2
        self.maximum = min(8, self.qtd)
        self.chosens = np.array([])
        self.ts = np.array([])
        self.costs = np.array([])

    def paths(self, chosen):
        S = np.where(chosen == 1)[1]
        # print("S = {}".format(S))
        perms = np.array(list(permutations(S)))

        K = np.zeros((len(perms), self.qtd))
        # print(perms)
        for i in range(len(perms)):
            K[i][perms[i][0]] = self.est[perms[i][0]]
            soma = self.est[perms[i][0]]
            # print(perms[i])
            for j in range(1, len(perms[i])):
                # print('pos', perms[i][j], 'soma: ', soma, 'anterior:', perms[i][j-1], 'tempo', self.p[perms[i][j-1]], 'est: ', self.est[perms[i][j]])
                soma = max(soma + self.p[perms[i][j-1]], self.est[perms[i][j]])
                K[i][perms[i][j]] = soma
            # print(K[i][:])
        # input(K)
        return K

    def valid_multipliers(self, path, t):
        if len(np.where(t > 0)[0]) < 2:  # not allowed clique with size < 2
            return False
        for i in range(len(path)):
            soma = 0
            for j in range(len(t)):
                soma += t[j] * path[i][j]
            if soma < 1:
                return False
        return True

    def mip(self, K, S):
        m = self.m
        m.clear()
        m.verbose = 0
        m.max_seconds = 90
        initial = time.time()
        t = [0] * self.qtd
        for i in S:
            t[i] = m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i))
        # t = np.array([m.add_var(var_type=CONTINUOUS, lb=0, name='t({})'.format(i)) for i in S])
        # c = m.add_var(var_type=CONTINUOUS, lb=0, name='c')
        m.objective = xsum(self.x[i] * t[i] for i in S)
        # print(range(len(K)))
        # print(range(len(K[0])))
        # print(t[3])
        # print(S)
        # print()
        # if len(K) > 2:
            # input(len(K))
        for i in range(len(K)):
            m += xsum(K[i][S[j]] * t[S[j]] for j in range(len(S))) >= 1, 'K({})'.format(i)
        # end = time.time()
        # print("MIP creation: {:4.3g}".format(end - initial))
        # initial = time.time()
        m.optimize()
        # end = time.time()
        # print("MIP solve: {:4.3g}".format(end - initial))
        # m.write('teste.lp')

        ts = np.zeros(self.qtd)
        if m.status != OptimizationStatus.OPTIMAL:
            return ts, 999999

        for i in S:
            ts[i] = t[i].x

        return ts

    def resolve_exact(self):
        qtd = np.math.factorial(self.qtd)/(np.math.factorial(self.minimum)*np.math.factorial(self.qtd - self.minimum))
        while qtd < self.max_exact and self.minimum <= self.maximum:
            comb = list(combinations(range(self.qtd), self.minimum))
            # print(comb)
            for c in comb:
                # input(c)
                chosen = np.array([np.zeros(self.qtd)])
                for i in c:
                    chosen[0][i] = 1
                # print(chosen)
                path = self.paths(chosen)
                t = self.mip(path, c)
                value = self.cost_function(t, chosen, path)
                self.exact += 1

                if round(value, 6) < 100:
                    # print(c, value)
                    self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
                    self.ts = np.vstack([self.ts, t]) if self.ts.size else t
                    self.costs = np.vstack([self.costs, value]) if self.costs.size else value
            self.minimum += 1
            if self.minimum < self.maximum:
                qtd = np.math.factorial(self.qtd) / (
                        np.math.factorial(self.minimum) * np.math.factorial(self.qtd - self.minimum))
            # print(self.minimum, self.maximum, qtd)
        # print('minimum: ', self.minimum)
        # print('qtd: ', qtd, 'max', self.max_exact)
        # print(self.chosens)
        # print(self.ts)
        # print(self.costs)
        # input()

    def random_init_mip(self):
        chosen = np.array([np.random.randint(2, size=self.qtd)])
        S = np.where(chosen == 1)[1]
        size_S = len(S)
        # key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
        # input(key)
        if size_S < self.minimum:
            while size_S < self.minimum:
                pos = np.random.randint(0, self.qtd)
                if chosen[0][pos] == 0:
                    chosen[0][pos] = 1
                    size_S += 1
            S = np.where(chosen == 1)[1]
        elif size_S > self.maximum:
            while size_S > self.maximum:
                pos = np.random.randint(0, size_S)
                if chosen[0][S[pos]] == 1:
                    chosen[0][S[pos]] = 0
                    size_S -= 1
            S = np.where(chosen == 1)[1]

        key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
        # input(key)
        path = self.paths(chosen)
        t = self.mip(path, S)
        value = self.cost_function(t, chosen, path)

        self.accepted[key] = (t, value)
        return chosen, t, value

    def random_neighbour_mip(self, chosen):
        neighborhood = np.random.choice([0, 1, 2, 3], 1, p=[0.25, 0.25, 0.25, 0.25])[0]
        # print(neighborhood)
        new_chosen = chosen.copy()
        S = np.where(chosen == 1)[1]
        # print(self.accepted)
        if neighborhood == 0:  # change some chosen at random
            # initial = time.time()
            # print('neighborhood', 0)
            pos = np.random.randint(0, self.qtd)
            new_chosen[0][pos] = 1 - chosen[0][pos]
            # print('pos', pos, 'new_chosen', new_chosen)
            # if is in limits
            S = np.where(new_chosen == 1)[1]
            # print('S', S)
            if len(S) < self.minimum or len(S) > self.maximum:
                key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
                # print('limit', key)
                # end = time.time()
                # print("time neighbor 0 max-minimum: {:4.3g}".format(end - initial))
                return chosen, self.accepted[key][0], self.accepted[key][1]

            key = '{}'.format(''.join(str(new_chosen[0][i]) for i in range(self.qtd)))
            # print(key)
            # if new configuration exists
            if key in self.accepted:
                # print('exists', new_chosen)
                # print(self.accepted[key])
                # end = time.time()
                # print("time neighbor 0 exists: {:4.3g}".format(end - initial))
                return new_chosen, self.accepted[key][0], self.accepted[key][1]
            # end = time.time()
            # print("time neighbor 0: {:4.3g}".format(end-initial))

        elif neighborhood == 1:  # change two chosens at random
            # print('neighborhood', 1)
            # initial = time.time()
            pos1 = np.random.randint(0, self.qtd)
            pos2 = np.random.randint(0, self.qtd)

            if pos1 == pos2:
                if pos1 == 0:
                    pos1 += 1
                elif pos1 == self.qtd:
                    pos1 -= 1

            new_chosen[0][pos1] = 1 - chosen[0][pos1]
            new_chosen[0][pos2] = 1 - chosen[0][pos2]
            # print('pos1', pos1, 'pos2', pos2, 'new_chosen', new_chosen)
            S = np.where(new_chosen == 1)[1]

            # if is in limits
            if len(S) < self.minimum or len(S) > self.maximum:
                key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
                # end = time.time()
                # print("time neighbor 1 maximum minimum: {:4.3g}".format(end - initial))
                return chosen, self.accepted[key][0], self.accepted[key][1]

            key = '{}'.format(''.join(str(new_chosen[0][i]) for i in range(self.qtd)))
            # if new configuration exists
            if key in self.accepted:
                # print('exists', key)
                # end = time.time()
                # print("time neighbor 1 exists: {:4.3g}".format(end - initial))
                return new_chosen, self.accepted[key][0], self.accepted[key][1]
            # print('dont existt', key)
            # end = time.time()
            # print("time neighbor 1: {:4.3g}".format(end-initial))

        elif neighborhood == 2:  # include one that is not in the solution if maximum not reached
            # initial = time.time()
            if len(S) == self.maximum:
                key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
                end = time.time()
                # print("time neighbor 2 maximum: {:4.3g}".format(end - initial))
                return chosen, self.accepted[key][0], self.accepted[key][1]

            S_not = np.where(chosen == 0)[1]
            # print('S_not', S_not)
            pos = S_not[np.random.randint(0, len(S_not))]
            # print('pos', pos)
            new_chosen[0][pos] = 1
            key = '{}'.format(''.join(str(new_chosen[0][i]) for i in range(self.qtd)))

            # if new configuration exists
            if key in self.accepted:
                # print('exists', key)
                # end = time.time()
                # print("time neighbor 2 exists: {:4.3g}".format(end - initial))
                return new_chosen, self.accepted[key][0], self.accepted[key][1]
            # print('dont existt', key)
            # end = time.time()
            # print("time neighbor 2: {:4.3g}".format(end-initial))

        else:  # remove one that is in solution
            # initial = time.time()
            if len(S) == self.minimum:
                key = '{}'.format(''.join(str(chosen[0][i]) for i in range(self.qtd)))
                # print('minimum', key)
                # end = time.time()
                # print("time neighbor 3 minimum : {:4.3g}".format(end - initial))
                return chosen, self.accepted[key][0], self.accepted[key][1]

            pos = S[np.random.randint(0, len(S))]
            # print('pos', pos)
            new_chosen[0][pos] = 0
            key = '{}'.format(''.join(str(new_chosen[0][i]) for i in range(self.qtd)))

            # if new configuration exists
            if key in self.accepted:
                # print('exists', key)
                # end = time.time()
                # print("time neighbor 3 exists: {:4.3g}".format(end - initial))
                return new_chosen, self.accepted[key][0], self.accepted[key][1]
            # end = time.time()
            # print("time neighbor 3: {:4.3g}".format(end-initial))
            # print('dont existt', key)
        # initial = time.time()
        S = np.where(new_chosen == 1)[1]
        K = self.paths(new_chosen)
        # end = time.time()
        # print("time path: {:4.3g}".format(end - initial))
        # initial = time.time()
        new_t = self.mip(K, S)
        # end = time.time()
        # print("time MIP: {:4.3g}".format(end - initial))
        # initial = time.time()
        new_cost = self.cost_function(new_t, new_chosen, K)
        # end = time.time()
        # print("time cost: {:4.3g}".format(end - initial))
        self.accepted[key] = (new_t, new_cost)
        return new_chosen, new_t, new_cost

    def random_init(self):
        chosen = np.array([np.random.randint(2, size=self.qtd)])
        t = np.array([np.zeros(self.qtd)])
        for i in range(self.qtd):
            t[0][i] = chosen[0][i]*np.random.random()
        path = self.paths(chosen)
        # valid = self.valid_multipliers(path, t)
        # print(chosen, t, path, valid)
        # input('teste')
        return chosen, t, path

    def cost_function(self, t, chosen, K):
        fo = np.sum(self.x * t)
        sum_t = np.sum(t)

        penal_min_size = max(0, 100000 * (2 - len(np.where(chosen == 1)[1])))
        lhs_K = 1 - np.sum(t * K, 1)

        penal_left_1 = lhs_K[lhs_K > 0].sum()
        # for i in range(len(path)):
        #     soma = 0
        #     lhs = np.sum(np.multiply(self.x, t))
        #     for j in range(len(t)):
        #         soma += t[j] * path[i][j]
        #     penal_left_1 += max(0, 1 - soma)

        # print(t, chosen)
        # print(K)
        # print(lhs_K, penal_min_size, fo, penal_left_1, sum_t)
        # print(100*fo, 100000*penal_left_1, 100000*penal_min_size, - 1e-8*sum_t)
        # input
        return 100*fo + 100000*penal_min_size + 100000*penal_left_1 - 1e-8*sum_t

    def initial_temperature_mip(self, chosen, t, cost):
        sum_delta = 0
        for i in range(10):
            initial = time.time()
            new_chosen, new_t, new_cost = self.random_neighbour_mip(chosen)
            end = time.time()
            # print(end-initial, cost, new_cost)
            sum_delta += abs(cost - new_cost)
        if round(sum_delta / 10, 10) == 0:
            return 100

        return sum_delta / 10


    def initial_temperature(self, chosen, t, path, cost):
        sum_delta = 0
        for i in range(10):
            new_chosen, new_t, new_path = self.random_neighbour(chosen, t, path)
            new_cost = self.cost_function(new_t, new_chosen, new_path)
            # print(cost, new_cost, abs(cost - new_cost))
            sum_delta += abs(cost - new_cost)
        return sum_delta

    def alpha(self, fraction: float):
        """ Example of temperature dicreasing as the process goes on."""
        return np.maximum(0.0001, 1 - fraction)

    def random_neighbour(self, chosen, t, path, force=False):
        # choose a element with probability p
        valids = np.where(chosen == 1)[1]
        if len(valids) < 2 or force:
            neighborhood = 0
        else:
            neighborhood = np.random.choice([0, 1, 2], 1, p=[0.005, 0.005, 0.99])[0]
        # print(neighborhood)
        new_chosen = chosen.copy()
        new_t = t.copy()
        changePath = False
        if neighborhood == 0:  # change some chosen element
            pos = np.random.randint(0, self.qtd)
            new_chosen[0][pos] = 1 - chosen[0][pos]
            new_t[0][pos] = new_chosen[0][pos] * np.random.random()
            # print(neighborhood, pos, chosen[0][pos], new_chosen[0][pos], new_t[0][pos])
            new_path = self.paths(new_chosen)
            changePath = True
        elif neighborhood == 1:  # change some value of a valid t
            # print(valids)
            pos = valids[np.random.randint(0, len(valids))]
            new_t[0][pos] = np.random.random()
            # print(neighborhood, t[0][pos], new_t[0][pos])
            new_path = path.copy()
        else:  # increase/decrease by 1-10% the value of some value of t
            pos = valids[np.random.randint(0, len(valids))]
            pct = np.random.choice([1, -1], 1)[0] * (np.random.randint(1, 16)/100)
            new_t_pos = new_t[0][pos] + pct * new_t[0][pos]
            if new_t_pos > 1:
                new_t_pos = new_t[0][pos] - pct * new_t[0][pos]
            new_t[0][pos] = new_t_pos
            # print(neighborhood, t[0][pos], new_t[0][pos])
            new_path = path.copy()

        # if changePath:
        #     if len(np.where(new_chosen == 1)[0]) < 2:  # not allowed clique with size < 2
        #         new_valid = False
        #
        # new_valid = self.valid_multipliers(new_path, new_t)

        # print('chosen: ', chosen, 'new: ', new_chosen)
        # print('t: ', t, 'new: ', new_t)
        # print('path: ', path, 'new: ', new_path)
        # print('valid:', valid, 'new: ', new_valid)
        # input()
        return new_chosen, new_t, new_path

    def acceptance_probability(self, cost, new_cost, temperature):
        if abs(round(new_cost, 8) - round(cost, 8)) < 1e-8:
            if self.debug: print("    - Acceptance probabilty = 1e-8...", end='')
            return 1e-8
        if round(new_cost, 8) < round(cost, 8):
            if self.debug: print("    - Acceptance probabilty = 1 as new_cost = {} < cost = {}...".format(new_cost, cost), end='')
            return 1
        else:
            prob = np.exp(- (new_cost - cost) / temperature)
            if self.debug: print("    - Acceptance probabilty = {:.3g}...".format(prob), end='')
            return prob

    def annealing_mip(self):
        """ Optimize the black-box function 'cost_function' with the simulated annealing algorithm."""
        solutions = 0
        self.resolve_exact()
        # print(self.minimum, self.maximum, self.exact)
        if self.minimum >= self.maximum:
            return self.chosens, self.ts, self.costs, self.minimum, self.maximum, self.exact
        chosen, t, cost = self.random_init_mip()
        # o quanto melhora ou piora de 10 vizinhos a média disso é a temperatura inicial
        if round(cost, 6) < 100:
            self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
            self.ts = np.vstack([self.ts, t]) if self.ts.size else t
            self.costs = np.vstack([self.costs, cost]) if self.costs.size else cost
            solutions += 1
        maxsteps = self.maxsteps
        T0 = self.initial_temperature_mip(chosen, t, cost)
        T = T0
        # input(T)
        worst = -1
        removal = -1

        rejected = 0
        step = 0
        force = False
        while step < maxsteps:  # and T > 1e-11:   #for step in range(maxsteps):
            step += 1
            initial = time.time()
            new_chosen, new_t, new_cost = self.random_neighbour_mip(chosen)
            # end = time.time()
            # print("Time neighbor: {:>4.3g}".format(end-initial))
            # print('custos: ', cost, new_cost)
            if self.debug: print("Step #{:>2}/{:>2} : T = {:>4.3g}, chosen = {}, t = {}, cost = {:>4.3g}, "
                            "new_chosen = {}, new_t = {}, new_cost = {:>4.3g} ...".format(step, maxsteps, T, chosen, t, cost, new_chosen, new_t, new_cost), end='')
            # print("Step #{:>2}/{:>2} : T = {:>4.3g}, accepted = {}, cost = {}, new_cost = {}".format(step, maxsteps, T, accepted, cost, new_cost), end='')
            # input()
            # initial = time.time()
            if self.acceptance_probability(cost, new_cost, T) > np.random.random():
                chosen, t, cost = new_chosen, new_t, new_cost
                rejected = 0
                if round(cost, 6) < 100:
                    # if solutions < self.max_solutions:
                        # print()
                    # initial = time.time()
                    if not any(( chosen == x ).all() for x in self.chosens):
                        # if solutions == 0 or np.sum(np.all(np.isclose(ts, t), axis=1)) == 0:  # if t is unique in self.ts
                        self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
                        self.ts = np.vstack([self.ts, t]) if self.ts.size else t
                        self.costs = np.vstack([self.costs, cost]) if self.costs.size else cost
                        # if cost > worst:
                        #     worst = cost
                        #     removal = solutions
                        # solutions += 1
                            # else:
                            #     print(ts)
                            #     print(t)
                            #     input()
                            # input(self.chosens)
                end = time.time()
                # print("Time accept: {:>4.3g}".format(end - initial))
                if self.debug: print("  ==> Accept it! Time = {:>4.3g}".format((end - initial)))
                # print("  ==> Accept it!")
            else:
                # print()
                end = time.time()
                # print("Time reject: {:>4.3g}".format(end - initial))
                if self.debug: print("  ==> Reject it... Time = {:>4.3g}".format(end-initial))
                # if T <= 1e-10:
                #     rejected += 1
                # if rejected > 0.2*maxsteps:
                #     # input()
                #     return self.chosens, self.ts, self.costs
            # initial = time.time()
            fraction = float(step / float(maxsteps))
            T = T * self.alpha(fraction)  # self.temperature(fraction)
            # end = time.time()
            # print("Temperature: {:>4.3g}. Time temperature: {:>4.3g}".format(T, end-initial))

            if T < 1e-10:
                return self.chosens, self.ts, self.costs, self.minimum, self.maximum, self.exact

            # input()
        # if self.debug: print('self.chosens', self.chosens)
        # input()
        return self.chosens, self.ts, self.costs, self.minimum, self.maximum, self.exact


    def annealing(self):
        """ Optimize the black-box function 'cost_function' with the simulated annealing algorithm."""
        solutions = 0
        chosen, t, path = self.random_init()
        cost = self.cost_function(t, chosen, path)
        self.chosens = np.array([])
        self.ts = np.array([])
        self.costs = np.array([])
        worst = -1
        removal = -1
        accepted = 0
        # o quanto melhora ou piora de 10 vizinhos a média disso é a temperatura inicial
        if cost < 100:
            self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
            self.ts = np.vstack([self.ts, t]) if self.ts.size else t
            self.costs = np.vstack([self.costs, cost]) if self.costs.size else cost
            solutions += 1
            accepted += 1
        maxsteps = self.maxsteps
        T0 = self.initial_temperature(chosen, t, path, cost)
        T = T0
        # input(T)
        rejected = 0
        step = 0
        force = False
        while step < maxsteps:  # and T > 1e-11:   #for step in range(maxsteps):
            step += 1
            new_chosen, new_t, new_path = self.random_neighbour(chosen, t, path, force)
            force = False
            new_cost = self.cost_function(new_t, new_chosen, new_path)
            # print('custos: ', cost, new_cost)
            if self.debug: print("Step #{:>2}/{:>2} : T = {:>4.3g}, chosen = {}, t = {}, cost = {:>4.3g}, "
                            "new_chosen = {}, new_t = {}, new_cost = {:>4.3g} ...".format(step, maxsteps, T, chosen, t, cost, new_chosen, new_t, new_cost), end='')
            # print("Step #{:>2}/{:>2} : T = {:>4.3g}, accepted = {}, cost = {}, new_cost = {}".format(step, maxsteps, T, accepted, cost, new_cost), end='')
            # input()
            if self.acceptance_probability(cost, new_cost, T) > np.random.random():
                chosen, t, cost, path = new_chosen, new_t, new_cost, new_path
                accepted += 1
                rejected = 0
                if cost < 100:
                    if solutions < self.max_solutions:
                        # print()
                        # if solutions == 0 or np.sum(np.all(np.isclose(ts, t), axis=1)) == 0:  # if t is unique in self.ts
                        self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
                        self.ts = np.vstack([self.ts, t]) if self.ts.size else t
                        self.costs = np.vstack([self.costs, cost]) if self.costs.size else cost
                        if cost > worst:
                            worst = cost
                            removal = solutions
                        solutions += 1
                        # else:
                        #     print(ts)
                        #     print(t)
                        #     input()
                    else:
                        if cost < worst: # checar se os multiplicadores são iguais também
                            # if np.sum(np.all(np.isclose(ts, t),
                            #                  axis=1)) == 0:  # if t is unique in self.ts
                            self.chosens[removal] = chosen
                            self.ts[removal] = t
                            self.costs[removal] = cost
                            worst = cost
                            #find next worst
                            for i in range(self.max_solutions):
                                if self.costs[i] > worst:
                                    removal = i
                                    worst = self.costs[i]
                if self.debug: print("  ==> Accept it!")
                # print("  ==> Accept it!")
            else:
                # print()
                if self.debug: print("  ==> Reject it...{}".format(rejected))
                if T <= 1e-10:
                    rejected += 1
                # if rejected > 0.2*maxsteps:
                #     # input()
                #     return self.chosens, self.ts, self.costs

            fraction = float(step / float(maxsteps))
            if T < 1e-10:
                T = T0
            else:
                T = T * self.alpha(fraction)  # self.temperature(fraction)

            # input()
        # if self.debug: print('self.chosens', self.chosens)
        # input()
        return self.chosens, self.ts, self.costs

    def LAHC(self, l):
        chosen, t, path = self.random_init()
        cost = self.cost_function(t, chosen, path)
        best_chosen, best_t, best_cost = chosen, t, cost
        self.chosens = np.array([])
        self.ts = np.array([])
        self.costs = np.array([])
        for i in range(l):
            self.chosens = np.vstack([self.chosens, chosen]) if self.chosens.size else chosen
            self.ts = np.vstack([ts, t]) if self.ts.size else t
            self.costs = np.vstack([costs, cost]) if self.costs.size else cost
        v = 0
        # print(self.chosens, self.ts, self.costs)
        for step in range(self.maxsteps):
            new_chosen, new_t, new_path = self.random_neighbour(chosen, t, path)
            new_cost = self.cost_function(new_t, new_chosen, new_path)
            if new_cost < best_cost or new_cost < self.costs[v]:
                self.costs[v] = new_cost
                self.chosens[v] = new_chosen
                self.ts[v] = new_t
                if new_cost < best_cost:
                    best_chosen, best_t, best_cost = new_chosen, new_t, new_cost
            v = (v + 1) % l
        # print(self.chosens, self.ts, self.costs)
        return best_chosen, best_t, best_cost

if __name__ == "__main__":
    x = [0, 0, 19, 10]
    p = [10, 1, 1, 9]
    est = [0, 0, 19, 10]
    clique = Clique(x, p, est, 10000)
    # clique.annealing_mip()
    # chosen, t, cost = clique.LAHC(10)
    # print('Chosen = {} ; t = {} ; cost = {}'.format(chosen, t, cost))
    # input()
    escolhidos, indices, custos, minimum, maximum, exact = clique.annealing_mip()
    print('minimum', minimum, 'maximum', maximum, 'exact', exact)
    for i in range(len(escolhidos)):
        print('Chosen = {} ; t = {} ; cost = {}'.format(escolhidos[i], indices[i], custos[i]))
