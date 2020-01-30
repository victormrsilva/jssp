import numpy as np
from itertools import permutations

class Clique:
    def __init__(self, x, p, est, maxsteps):
        self.qtd = len(x)
        self.p = p
        self.est = est
        self.x = x
        self.maxsteps = maxsteps
        self.debug = True
        self.max_solutions = 5

    def paths(self, chosen):
        S = np.where(chosen == 1)[1]
        perms = np.asarray(list(permutations(S)))

        K = np.zeros((len(perms), self.qtd))
        # print(perms)
        for i in range(len(perms)):
            soma = 0
            # print(perms[i])
            for j in range(len(perms[i])):
                if j == 0:
                    soma = self.est[perms[i][j]]
                    # print('pos', perms[i][j], 'soma: ', soma, 'est: ', self.est[perms[i][j]])
                else:
                    # print('pos', perms[i][j], 'soma: ', soma, 'anterior:', perms[i][j-1], 'tempo', self.p[perms[i][j-1]], 'est: ', self.est[perms[i][j]])
                    soma = max(soma + self.p[perms[i][j-1]], self.est[perms[i][j]])
                for aux in range(j, len(perms[j])):
                    K[i][perms[i][aux]] = soma
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
        # input()
        return 100*fo + 100000*penal_min_size + 100000*penal_left_1 - 1e-8*sum_t

    def temperature(self, fraction: float):
        """ Example of temperature dicreasing as the process goes on."""
        return 1000000*np.maximum(0.01, 1 - fraction)

    def random_neighbour(self, chosen, t, path):
        # choose a element with probability p
        valids = np.where(chosen == 1)[1]
        if len(valids) < 2:
            neighborhood = 0
        else:
            neighborhood = np.random.choice([0, 1, 2], 1, p=[0.02, 0.08, 0.9])[0]
        # print(neighborhood)
        new_chosen = chosen.copy()
        new_t = t.copy()
        changePath = False
        if neighborhood == 0:  # change some chosen element
            pos = np.random.randint(0, self.qtd)
            new_chosen[0][pos] = 1 - chosen[0][pos]
            new_t[0][pos] = new_chosen[0][pos] * np.random.random()
            # print(pos, chosen[pos], new_chosen[pos], new_t[pos])
            new_path = self.paths(new_chosen)
            changePath = True
        elif neighborhood == 1:  # change some value of a valid t
            # print(valids)
            pos = valids[np.random.randint(0, len(valids))]
            new_t[0][pos] = np.random.random()
            new_path = path.copy()
        else:  # increase/decrease by 1-10% the value of some value of t
            pos = valids[np.random.randint(0, len(valids))]
            new_t[0][pos] = min(new_t[0][pos] + (np.random.choice([1, -1], 1)[0] * (np.random.randint(1, 11)/100) * new_t[0][pos]), 1)
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

        if new_cost < cost:
            if self.debug: print("    - Acceptance probabilty = 1 as new_cost = {} < cost = {}...".format(new_cost, cost), end='')
            return 1
        else:
            prob = np.exp(- (new_cost - cost) / temperature)
            if self.debug: print("    - Acceptance probabilty = {:.3g}...".format(prob), end='')
            return prob

    def annealing(self):
        """ Optimize the black-box function 'cost_function' with the simulated annealing algorithm."""
        solutions = 0
        chosen, t, path = self.random_init()
        cost = self.cost_function(t, chosen, path)
        chosens = np.array([])
        ts = np.array([])
        costs = np.array([])
        worst = -1
        removal = -1
        if cost < 100:
            chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
            ts = np.vstack([ts, t]) if ts.size else t
            costs = np.vstack([costs, cost]) if costs.size else cost
            solutions += 1
        maxsteps = self.maxsteps
        for step in range(maxsteps):
            fraction = float(step / float(maxsteps))
            T = maxsteps * (1 - fraction) # self.temperature(fraction)
            new_chosen, new_t, new_path = self.random_neighbour(chosen, t, path)
            new_cost = self.cost_function(new_t, new_chosen, new_path)
            # print('custos: ', cost, new_cost)
            if self.debug: print("Step #{:>2}/{:>2} : T = {:>4.3g}, chosen = {}, t = {}, cost = {:>4.3g}, "
                            "new_chosen = {}, new_t = {}, new_cost = {:>4.3g} ...".format(step, maxsteps, T, chosen, t, cost, new_chosen, new_t, new_cost), end='')
            # input()
            if self.acceptance_probability(cost, new_cost, T) > np.random.random():
                chosen, t, cost, path = new_chosen, new_t, new_cost, new_path
                if cost < 100:
                    if solutions < self.max_solutions:
                        chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
                        ts = np.vstack([ts, t]) if ts.size else t
                        costs = np.vstack([costs, cost]) if costs.size else cost
                        if cost > worst:
                            worst = cost
                            removal = solutions
                        solutions += 1
                    else:
                        print(chosens)
                        print(costs)
                        print(ts)
                        print(chosens[removal], chosen)
                        print(ts[removal], t)
                        print(costs[removal], cost)
                        if cost < worst:
                            chosens[removal] = chosen
                            ts[removal] = t
                            costs[removal] = cost
                            worst = cost
                            #find next worst
                            for i in range(self.max_solutions):
                                if costs[i] > worst:
                                    removal = i
                                    worst = costs[i]
                            print(chosens)
                            print(costs)
                            print(ts)
                            print('new worst: ', worst, removal)
                            input()

                    # insert = True
                    # if chosens.ndim == 1:
                    #     if np.alltrue(chosens == new_chosen):  # if chosens are already inserted
                    #         # print(path)
                    #         if np.alltrue(new_t >= ts):  # if the chosen[i] is dominated
                    #             chosens = new_chosen
                    #             ts = new_t
                    #             costs = new_cost
                    #         insert = False
                    # else:
                    #     for i in range(len(chosens)):
                    #         if np.alltrue(chosens[i] == new_chosen):  # if chosens are already inserted
                    #             # print(path)
                    #             if np.alltrue(new_t >= ts[i]):  # if the chosen[i] is dominated
                    #                 chosens[i] = new_chosen
                    #                 ts[i] = new_t
                    #                 costs[i] = new_cost
                    #             insert = False
                    #             break
                    # if insert:
                    #     # print(path)
                    #     chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
                    #     ts = np.vstack([ts, t]) if ts.size else t
                    #     costs = np.vstack([costs, cost]) if costs.size else cost
                        # print(chosens[0], chosen, np.alltrue(chosens[0] == new_chosen))
                        # input(chosens)
                    # print('chosen: ', chosen )
                    # print('t: ', t)
                    # print('path: ', path)
                    # print('valid:', valid)
                    # input()
                if self.debug: print("  ==> Accept it!")
            else:
                if self.debug: print("  ==> Reject it...")
            # input()
        # if self.debug: print('chosens', chosens)
        return chosens, ts, costs



if __name__ == "__main__":
    x = [0, 0, 10, 19]
    p = [10, 1, 1, 10]
    est = [0, 0, 19, 10]
    clique = Clique(x, p, est, 1000)
    escolhidos, indices, custos = clique.annealing()
    for i in range(len(escolhidos)):
        print('Chosen = {} ; t = {} ; cost = {}'.format(escolhidos[i], indices[i], custos[i]))
