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

    def paths(self, chosen):
        S = np.where(chosen == 1)[0]
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
        chosen = np.random.randint(2, size=self.qtd)
        t = np.zeros(self.qtd)
        for i in range(self.qtd):
            t[i] = chosen[i]*np.random.random()
        path = self.paths(chosen)
        valid = self.valid_multipliers(path, t)
        # print(chosen, t, path, valid)
        # input('teste')
        return chosen, t, path, valid

    def cost_function(self, t, valid):
        total = np.sum(np.multiply(self.x, t))
        if valid:
            return total
        else:
            return total + 1000

    def temperature(self, fraction: float):
        """ Example of temperature dicreasing as the process goes on."""
        return np.maximum(0.01, np.minimum(1, 1 - fraction))

    def random_neighbour(self, chosen, t, path, valid):
        # choose a element with probability p
        valids = np.where(chosen == 1)[0]
        if len(valids) == 0:
            neighborhood = 0
        else:
            neighborhood = np.random.choice([0, 1, 2], 1, p=[0.2, 0.2, 0.6])[0]
        # print(neighborhood)
        new_chosen = chosen.copy()
        new_t = t.copy()
        changePath = False
        if neighborhood == 0:  # change some chosen element
            pos = np.random.randint(0, self.qtd)
            new_chosen[pos] = 1 - chosen[pos]
            new_t[pos] = new_chosen[pos] * np.random.random()
            # print(pos, chosen[pos], new_chosen[pos], new_t[pos])
            new_path = self.paths(new_chosen)
            changePath = True
        elif neighborhood == 1:  # change some value of a valid t
            # print(valids)
            pos = valids[np.random.randint(0, len(valids))]
            new_t[pos] = np.random.random()
            new_path = path.copy()
        else:  # increase/decrease by 1-10% the value of some value of t
            pos = valids[np.random.randint(0, len(valids))]
            new_t[pos] = new_t[pos] + (np.random.choice([1, -1], 1)[0] * (np.random.randint(1, 11)/100) * new_t[pos])
            if new_t[pos] > 1:
                new_t[pos] = 1
            new_path = path.copy()

        if changePath:
            if len(np.where(new_chosen == 1)[0]) < 2:  # not allowed clique with size < 2
                new_valid = False

        new_valid = self.valid_multipliers(new_path, new_t)

        # print('chosen: ', chosen, 'new: ', new_chosen)
        # print('t: ', t, 'new: ', new_t)
        # print('path: ', path, 'new: ', new_path)
        # print('valid:', valid, 'new: ', new_valid)
        # input()
        return new_chosen, new_t, new_path, new_valid

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
        chosen, t, path, valid = self.random_init()
        cost = self.cost_function(t, valid)
        chosens = np.array([])
        ts = np.array([])
        costs = np.array([])
        if valid and cost < 1:
            chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
            ts = np.vstack([ts, t]) if ts.size else t
            costs = np.vstack([costs, cost]) if costs.size else cost
        maxsteps = self.maxsteps
        for step in range(maxsteps):
            fraction = float(step / float(maxsteps))
            T = self.temperature(fraction)
            new_chosen, new_t, new_path, new_valid = self.random_neighbour(chosen, t, path, valid)
            new_cost = self.cost_function(new_t, new_valid)
            # print('custos: ', cost, new_cost)
            if self.debug: print("Step #{:>2}/{:>2} : T = {:>4.3g}, chosen = {}, t = {}, cost = {:>4.3g}, "
                            "new_chosen = {}, new_t = {}, new_cost = {:>4.3g} ...".format(step, maxsteps, T, chosen, t, cost, new_chosen, new_t, new_cost), end='')

            if self.acceptance_probability(cost, new_cost, T) > np.random.random():
                chosen, t, cost, path, valid = new_chosen, new_t, new_cost, new_path, new_valid
                if cost < 1:
                    insert = True
                    if chosens.ndim == 1:
                        if np.alltrue(chosens == new_chosen):  # if chosens are already inserted
                            # print(path)
                            if np.alltrue(new_t >= ts):  # if the chosen[i] is dominated
                                chosens = new_chosen
                                ts = new_t
                                costs = new_cost
                            insert = False
                    else:
                        for i in range(len(chosens)):
                            if np.alltrue(chosens[i] == new_chosen):  # if chosens are already inserted
                                # print(path)
                                if np.alltrue(new_t >= ts[i]):  # if the chosen[i] is dominated
                                    chosens[i] = new_chosen
                                    ts[i] = new_t
                                    costs[i] = new_cost
                                insert = False
                                break
                    if insert:
                        # print(path)
                        chosens = np.vstack([chosens, chosen]) if chosens.size else chosen
                        ts = np.vstack([ts, t]) if ts.size else t
                        costs = np.vstack([costs, cost]) if costs.size else cost
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
        if self.debug: print('chosens', chosens)
        return chosens, ts, costs



if __name__ == "__main__":
    x = [0, 0, 10, 19]
    p = [10, 1, 1, 10]
    est = [0, 0, 19, 10]
    clique = Clique(x, p, est)
    escolhidos, indices, custos = clique.annealing()
    for i in range(len(escolhidos)):
        print('Chosen = {} ; t = {} ; cost = {}'.format(escolhidos[i], indices[i], custos[i]))
