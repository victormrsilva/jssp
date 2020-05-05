from itertools import permutations, combinations
import random
import sys
import numpy as np
from mip.model import Model, xsum, maximize, VarList, ConstrList
from mip.constants import INTEGER, BINARY, CONTINUOUS, OptimizationStatus
import re

from collections import OrderedDict

class Teste:
    def __init__(self):
        self.model = Model(solver_name="gurobi")
        self.m = Model(name='chvatal', solver_name="gurobi")
        self.x = {}
        self.eps = 0.00000001

    def chvatal_gomory_zero(self):
        m = self.m
        model = self.model
        m.clear()
        m.verbose = 0
        m.max_seconds = 30
        eps = 0.0001
        delta = 0.01

        indexes = []
        ub = []
        zero = []
        for i in range(len(model.vars)):
            if abs(model.vars[i].x) > 0.00000001:
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

        u = [m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta,
                       name='u({})'.format(constrs[i].name)) for i in range(len(constrs))]
        a = {model.vars[i].name: m.add_var(var_type=INTEGER, lb=-100000, name='a({})'.format(model.vars[i])) for i in indexes}
        f = {i: m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f({})'.format(model.vars[i])) for i in indexes}
        a0 = m.add_var(var_type=INTEGER, lb=-100000, name='a0')
        f0 = m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f0')
        # obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

        newconstraints = {model.vars[i].name: 0 for i in indexes}
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
                elif var.name in newconstraintszero:
                    newconstraintszero[var.name] += mult * u[i] * coeff
                else:
                    input('errro')

            newconstraints['b'] += (mult * -1.0 * c.expr.const - modif_ub) * u[i]
        # print('non-UB/zero')
        # for (key, value) in newconstraints.items():
        #     print(key, ':', value)
        # print('zero')
        # for (key, value) in newconstraintszero.items():
        #     print(key, ':', value)
        # input()

        m.objective = maximize((xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0))
                                # - xsum(0.0001 * u[i] for i in range(len(u))))
        # m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[i] == newconstraints[model.vars[i].name] - a[model.vars[i].name], 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0, 'violated'
        m.optimize()
        m.write('cg_1.lp')
        print('obj: ', m.objective_value)
        # input('cg.lp')

        restrictions = []
        # print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                teste = 0
                for i in indexes:
                    teste += model.vars[i] * round(a[model.vars[i].name].xi(k), 10)
                var = 0
                b = 0
                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    # print('c:', c, ';u:', u[j].xi(k), ';rs: ', c.expr.const, ';mult: ', mult, ';b:', mult * (-1.0 * c.expr.const) * u[j].xi(k))
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                    for (v, coeff) in c.expr.expr.items():
                        # print(v, ': ', round(mult * coeff * u[j].xi(k), 10), 'var_antes: ', var, '<=', b)
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                # print('var: ', var, ' b: ', b)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                # print('eq: ', restr, ';soma: ', soma, ';b: ', b)
                # print('teste: ', teste, ';soma: ', xsum(model.vars[i].x * round(a[model.vars[i].name].xi(k), 10) for i in indexes),
                #       ';a0: ', round(a0.xi(k), 10))
                # print('==========================')
                if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                    print('eq: ', restr, '<=', b, 'diff: ', soma - b)
                    restrictions.append(restr <= b)
                    # restrictions.append(teste <= a0.xi(k))
        return restrictions

    def chvatal_gomory_zero_obj_penalties(self):
        m = self.m
        model = self.model
        m.clear()
        m.verbose = 0
        m.max_seconds = 30
        eps = 0.0001
        delta = 0.01

        indexes = []
        ub = []
        zero = []
        for i in range(len(model.vars)):
            if abs(model.vars[i].x) > 0.00000001:
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

        u = [m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta,
                       name='u({})'.format(constrs[i].name)) for i in range(len(constrs))]
        a = {model.vars[i].name: m.add_var(var_type=INTEGER, lb=-100000, name='a({})'.format(model.vars[i])) for i in indexes}
        f = {i: m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f({})'.format(model.vars[i])) for i in indexes}
        a0 = m.add_var(var_type=INTEGER, lb=-100000, name='a0')
        f0 = m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f0')
        # obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

        newconstraints = {model.vars[i].name: 0 for i in indexes}
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
                elif var.name in newconstraintszero:
                    newconstraintszero[var.name] += mult * u[i] * coeff
                else:
                    input('errro')

            newconstraints['b'] += (mult * -1.0 * c.expr.const - modif_ub) * u[i]
        # print('non-UB/zero')
        # for (key, value) in newconstraints.items():
        #     print(key, ':', value)
        # print('zero')
        # for (key, value) in newconstraintszero.items():
        #     print(key, ':', value)


        m.objective = maximize((xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0)
                                - xsum(0.0001 * u[i] for i in range(len(u))))
        # m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[i] == newconstraints[model.vars[i].name] - a[model.vars[i].name], 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0, 'violated'
        m.optimize()
        m.write('cg_1.lp')
        print('obj: ', m.objective_value)
        # input('cg.lp')

        restrictions = []
        # print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                teste = 0
                for i in indexes:
                    teste += model.vars[i] * round(a[model.vars[i].name].xi(k), 10)
                var = 0
                b = 0
                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    # print('c:', c, ';u:', u[j].xi(k), ';rs: ', c.expr.const, ';mult: ', mult, ';b:', mult * (-1.0 * c.expr.const) * u[j].xi(k))
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                    for (v, coeff) in c.expr.expr.items():
                        # print(v, ': ', round(mult * coeff * u[j].xi(k), 10), 'var_antes: ', var, '<=', b)
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                # print('var: ', var, ' b: ', b)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                # print('eq: ', restr, ';soma: ', soma, ';b: ', b)
                # print('teste: ', teste, ';soma: ', xsum(model.vars[i].x * round(a[model.vars[i].name].xi(k), 10) for i in indexes),
                #       ';a0: ', round(a0.xi(k), 10))
                # print('==========================')
                if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                    print('eq: ', restr, '<=', b, 'diff: ', soma - b)
                    restrictions.append(restr <= b)
                    # restrictions.append(teste <= a0.xi(k))
        return restrictions

    def chvatal_gomory_ub(self):
        m = self.m
        model = self.model
        m.clear()
        m.verbose = 0
        m.max_seconds = 30
        eps = 0.0001
        delta = 0.01

        indexes = []
        ub = []
        zero = []
        for i in range(len(model.vars)):
            if abs(model.vars[i].x) < 0.00000001:
                zero.append(i)
            elif abs(model.vars[i].x - model.vars[i].ub) < 0.00000001:
                ub.append(i)
            else:
                indexes.append(i)

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
        # obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

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
        print('non-UB/zero')
        for (key, value) in newconstraints.items():
            print(key, ':', value)
        print('UB')
        for (key, value) in newconstraintsUB.items():
            print(key, ':', value)
        print('zero')
        for (key, value) in newconstraintszero.items():
            print(key, ':', value)
        # input('check')
        m.objective = maximize((xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0))
                               # - xsum(0.0001 * u[i] for i in range(len(u))))
        # m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[i] == newconstraints[model.vars[i].name] - a[model.vars[i].name], 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0, 'violated'
        m.optimize()
        m.write('cg_1.lp')
        print('obj: ', m.objective_value)
        # input('check cg.lp')

        restrictions = []
        # print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                teste = 0
                for i in indexes:
                    teste += model.vars[i] * round(a[model.vars[i].name].xi(k), 10)
                var = 0
                b = 0

                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    print('c:', c, ';u:', u[j].xi(k), ';rs: ', c.expr.const, ';mult: ', mult, ';b:', mult * (-1.0 * c.expr.const) * u[j].xi(k))
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                    for (v, coeff) in c.expr.expr.items():
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                # print('var: ', var, ' b: ', b)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                print('eq: ', restr, ';soma: ', soma, ';b: ', b)
                print('teste: ', teste, ';soma: ', xsum(model.vars[i].x * round(a[model.vars[i].name].xi(k), 10) for i in indexes),
                      ';a0: ', round(a0.xi(k), 10))
                print('==========================')
                # if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                # soma_teste = 0
                # for i in indexes:
                    # soma_teste += model.vars[i].x * round(a[model.vars[i].name].xi(k), 10)
                if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                    print('eq: ', restr, '<=', round(b, 10), 'diff: ', round(soma - b,10))
                    # restrictions.append(teste <= round(a0.xi(k), 10))
                    restrictions.append(restr <= b)
        return restrictions


    def chvatal_gomory_ub_obj_penalties(self):
        m = self.m
        model = self.model
        m.clear()
        m.verbose = 0
        m.max_seconds = 30
        eps = 0.0001
        delta = 0.01

        indexes = []
        ub = []
        zero = []
        for i in range(len(model.vars)):
            if abs(model.vars[i].x) < 0.00000001:
                zero.append(i)
            elif abs(model.vars[i].x - model.vars[i].ub) < 0.00000001:
                ub.append(i)
            else:
                indexes.append(i)

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
        # obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

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
        print('non-UB/zero')
        for (key, value) in newconstraints.items():
            print(key, ':', value)
        print('UB')
        for (key, value) in newconstraintsUB.items():
            print(key, ':', value)
        print('zero')
        for (key, value) in newconstraintszero.items():
            print(key, ':', value)
        # input('check')
        m.objective = maximize((xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0)
                               - xsum(0.0001 * u[i] for i in range(len(u))))
        # m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[i] == newconstraints[model.vars[i].name] - a[model.vars[i].name], 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0, 'violated'
        m.optimize()
        m.write('cg_1.lp')
        print('obj: ', m.objective_value)
        # input('check cg.lp')

        restrictions = []
        # print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                teste = 0
                for i in indexes:
                    teste += model.vars[i] * round(a[model.vars[i].name].xi(k), 10)
                var = 0
                b = 0

                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    print('c:', c, ';u:', u[j].xi(k), ';rs: ', c.expr.const, ';mult: ', mult, ';b:', mult * (-1.0 * c.expr.const) * u[j].xi(k))
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                    for (v, coeff) in c.expr.expr.items():
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                # print('var: ', var, ' b: ', b)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                print('eq: ', restr, ';soma: ', soma, ';b: ', b)
                print('teste: ', teste, ';soma: ', xsum(model.vars[i].x * round(a[model.vars[i].name].xi(k), 10) for i in indexes),
                      ';a0: ', round(a0.xi(k), 10))
                print('==========================')
                # if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                # soma_teste = 0
                # for i in indexes:
                    # soma_teste += model.vars[i].x * round(a[model.vars[i].name].xi(k), 10)
                if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                    print('eq: ', restr, '<=', round(b, 10), 'diff: ', round(soma - b, 10))
                    # restrictions.append(teste <= round(a0.xi(k), 10))
                    restrictions.append(restr <= b)
        return restrictions

    def chvatal_gomory_no_tune(self):
        m = self.m
        model = self.model
        model.verbose = 0
        m.clear()
        m.verbose = 0
        m.max_seconds = 30
        eps = 0.0001
        delta = 0.01

        indexes = []

        for i in range(len(model.vars)):
            indexes.append(i)

        for i in range(len(model.vars)):
            print(' {}: {}'.format(model.vars[i].name, round(model.vars[i].x, 10)), end='')
            if abs(model.vars[i].x - model.vars[i].ub) < eps:
                print('*', end='')

        print()
        constrs = [c for c in model.constrs if not c.name.startswith('gomory')]
        # for i in range(len(constrs)):
        #     print(i, constrs[i])

        u = [m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta,
                       name='u({})'.format(constrs[i].name)) for i in range(len(constrs))]
        a = {model.vars[i].name: m.add_var(var_type=INTEGER, lb=-100000, name='a({})'.format(model.vars[i].name)) for i in indexes}
        f = {model.vars[i].name: m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f({})'.format(model.vars[i].name)) for i in indexes}
        a0 = m.add_var(var_type=INTEGER, lb=-100000, name='a0')
        f0 = m.add_var(var_type=CONTINUOUS, lb=0, ub=1 - delta, name='f0')
        # obj = m.add_var(var_type=CONTINUOUS, lb=-100000, name='obj')

        newconstraints = {model.vars[i].name: 0 for i in indexes}
        newconstraints['b'] = 0
        # print(newconstraints)
        for i in range(len(constrs)):
            mult = 1
            c = model.constrs[i]
            if c.expr.sense == '>':
                mult = -1
            for (var, coeff) in c.expr.expr.items():
                newconstraints[var.name] += mult * u[i] * coeff
            newconstraints['b'] += (mult * -1.0 * c.expr.const) * u[i]

        m.objective = maximize(xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0)
        # m += obj == xsum(model.vars[i].x * a[model.vars[i].name] for i in indexes) - a0, 'objective'
        for i in indexes:
            m += f[model.vars[i].name] == newconstraints[model.vars[i].name] - a[model.vars[i].name], \
                 'f_{}'.format(model.vars[i].name)
        m += f0 == newconstraints['b'] - a0, 'f'
        # m += xsum(a[model.vars[i].name] * model.vars[i].x for i in indexes) >= a0 + 0.001, 'violated'
        # m.write('cg_2.lp')
        m.optimize()
        m.write('cg_4.lp')
        print('obj: ', m.objective_value)
        # input('cg.lp')

        restrictions = []
        print(m.num_solutions, m.status)
        if m.status == OptimizationStatus.OPTIMAL or m.status == OptimizationStatus.FEASIBLE:
            for k in range(m.num_solutions):
                teste = 0
                for i in indexes:
                    # print(model.vars[i].name, round(a[model.vars[i].name].xi(k), 10))
                    teste += model.vars[i] * round(a[model.vars[i].name].xi(k), 10)
                var = 0
                b = 0
                for j in range(len(constrs)):
                    mult = 1
                    c = constrs[j]
                    if c.expr.sense == '>':
                        mult = -1
                    # print('c:', c, ';u:', u[j].xi(k), ';rs: ', c.expr.const, ';mult: ', mult, ';b:', mult * (-1.0 * c.expr.const) * u[j].xi(k))
                    b += mult * (-1.0 * c.expr.const) * u[j].xi(k)
                    for (v, coeff) in c.expr.expr.items():
                        # print(v, ': ', round(mult * coeff * u[j].xi(k), 10), 'var_antes: ', var, '<=', b)
                        var += round(mult * coeff * u[j].xi(k), 10) * v
                b = np.floor(round(b, 10))
                restr = 0
                soma = 0
                # print('var: ', var, ' b: ', b)
                for (v, coeff) in var.expr.items():
                    restr += np.floor(round(coeff, 10)) * v
                    soma += np.floor(round(coeff, 10)) * v.x

                # print('eq: ', restr, 'soma; ', soma, 'b: ', b, 'soma > b?', soma > b)
                # print('teste: ', teste, ';soma:', xsum(model.vars[i].x * round(a[model.vars[i].name].xi(k), 10) for i in indexes),
                #       ';a0:', a0.xi(k))
                # print('==========================')
                if round(soma, 10) > round(b, 10) + 0.0001 and len(restr.expr.items()) > 0:
                    print('eq: ', restr, '<=', b, 'diff: ', soma - b)
                    restrictions.append(restr <= b)
        return restrictions


    def teste(self, lp: str, funcao: int):
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
        obj_inicial = self.model.objective_value
        print('{} - obj: {} ; atual = {}'.format(lp, obj, obj_inicial))

        if funcao == 1:
            cg = self.chvatal_gomory_zero()
        elif funcao == 2:
            cg = self.chvatal_gomory_zero_obj_penalties()
        elif funcao == 3:
            cg = self.chvatal_gomory_ub()
        elif funcao == 4:
            cg = self.chvatal_gomory_ub_obj_penalties()
        else:
            cg = self.chvatal_gomory_no_tune()
        totalCG = 0
        print('{} - obj: {} ; atual = {} ; CG encontrado: {}'.format(lp, obj, self.model.objective_value, len(cg)))
        while len(cg) > 0:
            for r in cg:
                print(r)
                totalCG += 1
                self.model += r, 'gomory{}'.format(totalCG)
            self.model.optimize()
            self.model.write('teste.lp')
            # if obj < self.model.objective_value:
            if funcao == 1:
                cg = self.chvatal_gomory_zero()
            elif funcao == 2:
                cg = self.chvatal_gomory_zero_obj_penalties()
            elif funcao == 3:
                cg = self.chvatal_gomory_ub()
            elif funcao == 4:
                cg = self.chvatal_gomory_ub_obj_penalties()
            else:
                cg = self.chvatal_gomory_no_tune()
            print('{} - obj: {} ; atual = {} ; CG encontrado: {}'.format(lp, obj, self.model.objective_value, len(cg)))

        self.model.optimize()
        # input('{} - otimo = {} ; inicial = {} ; atual = {} ; TotalCG: {}'.format(lp, obj, obj_inicial,
        #                                                                          self.model.objective_value, totalCG))
        return obj, obj_inicial, self.model.objective_value, totalCG
        # if obj < self.model.objective_value:
            # input('{} - obj: {} ; atual = {} ; CG encontrado: {}'.format(lp, obj, self.model.objective_value, len(cg)))

    def criaTestes(self):
        maxlp = 50
        # self.model.verbose = 0
        lp = 0
        while lp < maxlp:
            # print('/***************/')
            self.model.clear()
            self.model.verbose = 0
            self.x.clear()
            numVars = random.randint(2, 5)
            numConstr = random.randint(2, 5)
            for i in range(numVars):
                self.x['x{}']=self.model.add_var('x{}'.format(i), lb=-random.randint(0,5), ub=random.randint(1,5), var_type=INTEGER)

            obj = 0
            for j in range(numVars):
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
                        lin += rand * coef * self.x['x{}'.format(j)]
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
                print(self.model.objective_value, obj)
                if self.model.objective_value < obj:
                    cg = self.chvatal_gomory_zero_obj_penalties()
                    totalCG = 0
                    while len(cg) > 0:
                        for r in cg:
                            print(r)
                            totalCG += 1
                            self.model += r, 'gomory{}'.format(totalCG)
                        self.model.optimize()
                        self.model.write('teste.lp')
                        if self.model.status != OptimizationStatus.OPTIMAL:
                            input('erro. checar teste.lp')
                        if obj < self.model.objective_value:
                            print(obj, self.model.objective_value)
                            input('valor da fo invalido. checar teste.lp')
                        cg = self.chvatal_gomory_zero_obj_penalties()
                    # if not self.integerSol():
                    if totalCG > 0:
                        print(totalCG)
                        print('teste{}.lp'.format(lp))
                        lp = lp + 1
                        if round(obj) - round(self.model.objective_value) < -0.000001:
                            print(obj, ' ', self.model.objective_value)
                            input('veja lp {}'.format(lp - 1))
                # else:
                #     input('veja')
            # print('==============')

    def criaTestes_2vars(self):
        maxlp = 50
        # self.model.verbose = 0
        lp = 0
        while lp < maxlp:
            # print('/***************/')
            self.model.clear()
            self.model.verbose = 0
            self.x.clear()
            positive = []
            negative = []
            numVars = random.randint(2, 5)
            numConstr = random.randint(2, 5)
            for i in range(numVars):
                while True:
                    lb = random.randint(0, 5)
                    ub = random.randint(0, 5)
                    if lb != 0 or ub != 0:
                        break

                if lb == 0:
                    self.x['x{}'.format(i)] = self.model.add_var('x{}'.format(i), lb=0, ub=ub, var_type=INTEGER)
                    positive.append(i)
                else:
                    self.x['x{}p'.format(i)] = self.model.add_var('x{}p'.format(i), lb=0, ub=ub, var_type=INTEGER)
                    self.x['x{}n'.format(i)] = self.model.add_var('x{}n'.format(i), lb=0, ub=lb, var_type=INTEGER)
                    negative.append(i)

            obj = 0
            for j in range(numVars):
                rand = random.randint(0, 1)
                if rand == 0:
                    rand = -1
                coef = random.randint(0, 5)
                if j in positive:
                    obj += rand * coef * self.x['x{}'.format(j)]
                else:
                    obj += rand * coef * (self.x['x{}p'.format(j)] - self.x['x{}n'.format(j)])


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
                        if j in positive:
                            lin += rand * coef * self.x['x{}'.format(j)]
                        else:
                            lin += rand * coef * (self.x['x{}p'.format(j)] - self.x['x{}n'.format(j)])
                # print(lin)
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
            # self.model.write('teste.lp'.format(lp))
            # input('teste.lp')
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
                print(self.model.objective_value, obj)
                if self.model.objective_value < obj:
                    cg = self.chvatal_gomory_no_tune()
                    totalCG = 0
                    while len(cg) > 0:
                        for r in cg:
                            print(r)
                            totalCG += 1
                            self.model += r, 'gomory{}'.format(totalCG)
                        self.model.optimize()
                        self.model.write('teste.lp')
                        if self.model.status != OptimizationStatus.OPTIMAL:
                            input('erro. checar teste.lp')
                        if obj < self.model.objective_value:
                            print(obj, self.model.objective_value)
                            input('valor da fo invalido. checar teste.lp')
                        cg = self.chvatal_gomory_no_tune()
                    # if not self.integerSol():
                    if totalCG > 0:
                        print(totalCG)
                        print('teste{}.lp'.format(lp))
                        lp = lp + 1
                        if round(obj) - round(self.model.objective_value) < -0.000001:
                            print(obj, ' ', self.model.objective_value)
                            input('veja lp {}'.format(lp - 1))
                # else:
                #     input('veja')
            # print('==============')


    def integerSol(self):
        for v in self.model.vars:
            print(v, v.x)
            if not self.is_whole(v.x):
                return False
        return True

    def is_whole(self, f):
        return abs(f - round(f)) < abs(self.eps)

if __name__ == "__main__":
    # for i in range(17):
    #     t = Teste()
    #     t.teste('teste{}.lp'.format(i))
    #     del t
    #     input('/***************************************************/')
    cg = int(sys.argv[2])
    t = Teste()
    results = []
    # for i in range(50):
    #     results.append(t.teste('teste{}.lp'.format(i), cg))
    print('========================')
    results.append(t.teste(sys.argv[1], cg))
    for i in range(len(results)):
        print(results[i])
    # t.teste(sys.argv[1], cg)
    # t.teste('teste0.lp')
    # t.teste('rcpsp7.lp')
    # t.teste('vi_33_28_cuts.lp')
    # t.teste('vi4439_modelMCLin.lp')
    # t.teste('vi3328_modelMCLin.lp')
    # t.criaTestes_2vars()

