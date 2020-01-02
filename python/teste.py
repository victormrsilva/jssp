from itertools import permutations, combinations
import random

import numpy as np
from mip.model import Model, xsum, maximize, VarList, ConstrList
from mip.constants import INTEGER, BINARY, CONTINUOUS, OptimizationStatus

from collections import OrderedDict

class Teste:
    def __init__(self):
        self.model = Model()
        self.m = Model(name='chvatal', solver_name="gurobi")
        self.x = []
        self.eps = 0.00000001

    def chvatal_gomory(self):
        m = self.m
        model = self.model
        m.clear
        m.verbose = 1
        m.max_seconds = 60
        eps = 0.0001
        delta = 0.01

        indexes = []
        ub = []
        zero = []
        for i in range(len(model.vars)):
            if abs(model.vars[i].x - model.vars[i].ub) < 0.00000001:
                ub.append(i)
            elif abs(model.vars[i].x) > 0.00000001:
                indexes.append(i)
            else:
                zero.append(i)

        for i in range(len(model.vars)):
            print(model.vars[i].name, model.vars[i].x, end='')
            if abs(model.vars[i].x - model.vars[i].ub) < eps:
                print(' *')
            else:
                print()


        constrs = [c for c in model.constrs if not c.name.startswith('gomory')]
        # for i in range(len(constrs)):
        #     print(i, constrs[i])

        u = [m.add_var(var_type=CONTINUOUS, lb=(0, -1 + delta)[constrs[i].expr.sense == '='], ub=1 - delta,
                       name='u({})'.format(constrs[i].name)) for i in range(len(constrs))]
        a = {model.vars[i].name: m.add_var(var_type=INTEGER, lb=-100000, name='a({})'.format(model.vars[i])) for i in indexes}
        f = {i: m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f({})'.format(model.vars[i])) for i in indexes}
        a0 = m.add_var(var_type=INTEGER, lb=-100000, name='a0')
        f0 = m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f0')
        obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

        newconstraints = {model.vars[i].name: 0 for i in indexes}
        newconstraintsUB = {model.vars[i].name: 0 for i in ub}
        newconstraintszero = {model.vars[i].name: 0 for i in zero}
        newconstraints['b'] = 0

        for i in range(len(constrs)):
            mult = 1
            c = constrs[i]
            modif_ub = 0
            if c.expr.sense == '>':
                mult = -1
            for (var, coeff) in c.expr.expr.items():
                if var.name in newconstraints:
                    newconstraints[var.name] += mult * u[i] * coeff
                elif var.name in newconstraintsUB:
                    newconstraintsUB[var.name] += mult * u[i] * coeff * -1
                    # print(var.x, mult, coeff, -1)
                    modif_ub += var.x * mult * coeff
                    # print(modif_ub)
                elif var.name in newconstraintszero:
                    newconstraintszero[var.name] += mult * u[i] * coeff

            newconstraints['b'] += (mult * -1.0 * c.expr.const - modif_ub) * u[i]
        print('non-UB non-zero')
        for (key, value) in newconstraints.items():
            print(key, ':', value)
        print('UB')
        for (key, value) in newconstraintsUB.items():
            print(key, ':', value)
        print('zero')
        for (key, value) in newconstraintszero.items():
            print(key, ':', value)
        # input()

        m.objective = maximize(obj - xsum(0.0001 * u[i] for i in range(len(u)) if abs(constrs[i].slack) <= eps))
        m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[i] == newconstraints[model.vars[i].name] - a[model.vars[i].name], 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0, 'violated'
        m.write('cg.lp')
        # input('cg.lp')
        m.optimize()


        restrictions = []
        print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                # soma = xsum(a[vars[i].name].xi(k) * vars[i].x for i in indexes)
                # var = xsum(a[model.vars[i].name].xi(k) * model.vars[i] for i in indexes)
                # print(soma.const, a0.xi(k))
                # varUB = 0
                # for i in ub:
                #     coefVar = 0
                #     print(model.vars[i].name, ': ', newconstraintsUB[model.vars[i].name], '|',  end='')
                #     if not isinstance(newconstraintsUB[model.vars[i].name], int):
                #         for (var, coeff) in newconstraintsUB[model.vars[i].name].expr.items():
                #             if abs(var.xi(k)) > 0.00000001:
                #                 print(' ', var.name, ': ', round(var.xi(k), 10), ' (', coeff * round(var.xi(k), 10), ')', end='')
                #             coefVar += coeff * round(var.xi(k), 10)
                #     print()
                #     varUB += np.floor(coefVar) * model.vars[i]
                # varZero = 0
                # for i in zero:
                #     coefVar = 0
                #     print(model.vars[i].name, ': ', newconstraintszero[model.vars[i].name], '|', end='')
                #     if not isinstance(newconstraintszero[model.vars[i].name], int):
                #         for (var, coeff) in newconstraintszero[model.vars[i].name].expr.items():
                #             if abs(var.xi(k)) > 0.00000001:
                #                 print(' ', var.name, ': ', round(var.xi(k), 10), ' (', coeff * round(var.xi(k), 10), ')', end='')
                #             coefVar += coeff * round(var.xi(k), 10)
                #     print()
                #     varZero += np.floor(coefVar) * model.vars[i]
                # varNormal = 0
                # for i in indexes:
                #     coefVar = 0
                #     print(model.vars[i].name, ': ', newconstraints[model.vars[i].name], '|', end='')
                #     if not isinstance(newconstraints[model.vars[i].name], int):
                #         for (var, coeff) in newconstraints[model.vars[i].name].expr.items():
                #             if abs(var.xi(k)) > 0.00000001:
                #                 print(' ', var.name, ': ', round(var.xi(k), 10), ' (', coeff * round(var.xi(k), 10), ')', end='')
                #             coefVar += coeff * round(var.xi(k), 10)
                #         print()
                #     varNormal += np.round(coefVar) * model.vars[i]

                # var = xsum(round(a[model.vars[i].name].xi(k)) * model.vars[i] for i in indexes)
                # soma = 0
                # var = varNormal + varUB + varZero
                var = 0
                b = 0
                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    # print('c:', c, 'u:', u[j].xi(k), 'rs: ', c.expr.const, 'mult: ', mult)
                    for (v, coeff) in c.expr.expr.items():
                        # print(v, mult, u[j].xi(k))
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                print('var: ', var)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                print('eq: ', restr, 'soma; ', soma, 'b: ', b)
                print('==========================')
                if soma > b and len(restr.expr.items()) > 0:
                    restrictions.append(restr <= b)
                # soma = xsum(a[model.vars[i].name].xi(k) * model.vars[i].x for i in indexes)
                # print(var, round(soma, 8), a0.xi(k))
                # if round(soma, 8) > a0.xi(k) and len(var.expr.items()) > 0:
                #     restrictions.append(var <= a0.xi(k))
        return restrictions

    def teste(self, lp: str):
        self.model.clear()
        self.model.verbose = 0
        self.model.read(lp)
        self.model.optimize()
        obj = self.model.objective_value
        self.model.clear()
        self.model.verbose = 0
        self.model.read(lp)
        self.model.relax()
        self.model.optimize()

        cg = self.chvatal_gomory()
        totalCG = 0
        print('{} - obj: {} ; atual = {} '.format(lp, obj, self.model.objective_value))
        input(cg)
        while len(cg) > 0:
            for r in cg:
                print(r)
                totalCG += 1
                self.model += r, 'gomory{}'.format(totalCG)
            self.model.optimize()
            self.model.write('teste.lp')
            # if obj < self.model.objective_value:
            cg = self.chvatal_gomory()
            print('{} - obj: {} ; atual = {} '.format(lp, obj, self.model.objective_value))
            input(cg)

        self.model.optimize()
        if obj < self.model.objective_value:
            input('{} - obj: {} ; atual = {} '.format(lp, obj, self.model.objective_value))

    def criaTestes(self):
        maxlp = 50
        # self.model.verbose = 0
        lp = 0
        while lp < maxlp:
            # print('/***************/')
            self.model.clear()
            self.model.verbose = 0
            del self.x[:]
            numVars = random.randint(2, 5)
            numConstr = random.randint(2, 5)
            for i in range(numVars):
                self.x.append(self.model.add_var('x{}'.format(i), lb=0, ub=random.randint(1,5), var_type=INTEGER))

            obj = 0
            for j in range(numVars):
                rand = random.randint(0, 1)
                rand = random.randint(0, 1)
                if rand == 0:
                    rand = -1
                coef = random.randint(0, 5)
                obj += rand * coef * self.x[j]

            # print('obj: ', obj)
            self.model.objective = obj

            i = 0
            while i < numConstr:
                # print(i)
                i = i + 1
                lin = 0
                rand = 0
                for j in range(0, numVars):
                    rand = random.randint(0, 1)
                    if rand == 0:
                        rand = -1
                    coef = random.randint(0, 5)

                    if coef > 0:
                        lin += rand * coef * self.x[j]
                rand = random.randint(0, 2)
                if isinstance(lin, int) or len(lin.expr) < 2:
                    # print(lin)
                    i = i - 1
                else:
                    if rand == 0:
                        lin = lin == random.randint(0, 5 * numVars)
                    elif rand == 1:
                        lin = lin <= random.randint(0, 5 * numVars)
                    else:
                        lin = lin >= -1 * random.randint(0, 1) * random.randint(0, 5 * numVars)
                    # print('c', i, ': ', lin)
                    self.model += lin
            self.model.write('teste{}.lp'.format(lp))
            self.model.optimize()
            obj = 0
            # print(self.model.status)
            if self.model.status == OptimizationStatus.OPTIMAL:
                obj = self.model.objective_value
            else:
                continue
            self.model.relax()
            self.model.write('teste{}_relax.lp'.format(lp))
            self.model.optimize()
            # print(self.model.status)
            if self.model.status == OptimizationStatus.OPTIMAL:
                # print(self.model.objective_value, obj)
                if self.model.objective_value < obj:
                    cg = self.chvatal_gomory()
                    totalCG = 0
                    while len(cg) > 0:
                        for r in cg:
                            print(r)
                            totalCG += 1
                            self.model += r, 'gomory{}'.format(totalCG)
                        self.model.optimize()
                        self.model.write('teste.lp')
                        # if obj < self.model.objective_value:
                        cg = self.chvatal_gomory()
                    # if not self.integerSol():
                    if totalCG > 0:
                        print('teste{}.lp'.format(lp))
                        lp = lp + 1
                # else:
                    # input('veja')
            # print('==============')

    def integerSol(self):
        for v in self.model.vars:
            print(v, v.x)
            if not self.is_whole(v.x):
                return False
        return True

    def is_whole(self, f):
        return abs(f - round(f)) < abs(self.eps)

def test1():
    model = Model()

    x1 = model.add_var(var_type=BINARY, name='x(1)')
    x2 = model.add_var(var_type=INTEGER, name='x(2)', lb=0, ub=2)
    x3 = model.add_var(var_type=INTEGER, name='x(3)', lb=0, ub=2)
    model.objective = maximize(x1 + 2*x2 + 2*x3)
    model += 5*x1 + 5*x2 + x3 <= 11, 'c(1)'
    model += x1 + x2 + x3 >= 3, 'c(2)'

    model.relax()

    model.optimize()

    for var in model.vars:
        print(var.name, var.x)
    for c in model.constrs:
        print(c)

    cg = chvatal_gomory(model)
    totalCG = 0
    # input(cg)
    while len(cg) > 0:
        for r in cg:
            print(r)
            totalCG += 1
            model += r, 'gomory{}'.format(totalCG)
        model.optimize()
        model.write('teste.lp')
        input('pausa')
        cg = chvatal_gomory(model)


if __name__ == "__main__":
    # for i in range(50):
        # if i in [1, 3]:
        #     continue
        # t = Teste()
        # t.teste('teste{}.lp'.format(i))
        # del t
        # input('/***************************************************/')
    t = Teste()
    # t.teste('rcpsp6.lp')
    t.teste('vi_33_28_cuts.lp')
    # t.teste('vi4439_modelMCLin.lp')
    # t.teste('vi3328_modelMCLin.lp')
    # t.criaTestes()


