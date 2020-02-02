import numpy as np
from itertools import permutations

class Clique:
    def __init__(self, x, p, est, maxsteps):
        self.qtd = len(x)
        self.p = p
        self.est = est
        self.x = x
        self.maxsteps = maxsteps
        self.debug = False
        self.max_solutions = 50

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
        # input
        return 100*fo + 100000*penal_min_size + 100000*penal_left_1 - 1e-8*sum_t

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
        accepted = 0
        # o quanto melhora ou piora de 10 vizinhos a média disso é a temperatura inicial
        if cost < 100:
            chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
            ts = np.vstack([ts, t]) if ts.size else t
            costs = np.vstack([costs, cost]) if costs.size else cost
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
                        # if solutions == 0 or np.sum(np.all(np.isclose(ts, t), axis=1)) == 0:  # if t is unique in ts
                        chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
                        ts = np.vstack([ts, t]) if ts.size else t
                        costs = np.vstack([costs, cost]) if costs.size else cost
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
                            #                  axis=1)) == 0:  # if t is unique in ts
                            chosens[removal] = chosen
                            ts[removal] = t
                            costs[removal] = cost
                            worst = cost
                            #find next worst
                            for i in range(self.max_solutions):
                                if costs[i] > worst:
                                    removal = i
                                    worst = costs[i]
                if self.debug: print("  ==> Accept it!")
                # print("  ==> Accept it!")
            else:
                # print()
                if self.debug: print("  ==> Reject it...{}".format(rejected))
                if T <= 1e-10:
                    rejected += 1
                # if rejected > 0.2*maxsteps:
                #     # input()
                #     return chosens, ts, costs

            fraction = float(step / float(maxsteps))
            if T < 1e-10:
                T = T0
            else:
                T = T * self.alpha(fraction)  # self.temperature(fraction)

            # input()
        # if self.debug: print('chosens', chosens)
        # input()
        return chosens, ts, costs

    def LAHC(self, l):
        chosen, t, path = self.random_init()
        cost = self.cost_function(t, chosen, path)
        best_chosen, best_t, best_cost = chosen, t, cost
        chosens = np.array([])
        ts = np.array([])
        costs = np.array([])
        for i in range(l):
            chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
            ts = np.vstack([ts, t]) if ts.size else t
            costs = np.vstack([costs, cost]) if costs.size else cost
        v = 0
        # print(chosens, ts, costs)
        for step in range(self.maxsteps):
            new_chosen, new_t, new_path = self.random_neighbour(chosen, t, path)
            new_cost = self.cost_function(new_t, new_chosen, new_path)
            if new_cost < best_cost or new_cost < costs[v]:
                costs[v] = new_cost
                chosens[v] = new_chosen
                ts[v] = new_t
                if new_cost < best_cost:
                    best_chosen, best_t, best_cost = new_chosen, new_t, new_cost
            v = (v + 1) % l
        # print(chosens, ts, costs)
        return best_chosen, best_t, best_cost

if __name__ == "__main__":
    x = [0, 0, 10, 19]
    p = [10, 1, 1, 10]
    est = [0, 0, 19, 10]
    clique = Clique(x, p, est, 10000)
    chosen, t, cost = clique.LAHC(10)
    print('Chosen = {} ; t = {} ; cost = {}'.format(chosen, t, cost))
    input()
    escolhidos, indices, custos = clique.annealing()
    for i in range(len(escolhidos)):
        print('Chosen = {} ; t = {} ; cost = {}'.format(escolhidos[i], indices[i], custos[i]))
