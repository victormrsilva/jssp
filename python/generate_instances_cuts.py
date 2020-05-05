from sys import argv
import random
from shutil import copyfile

import time
from JSSPInstance import JSSPInstance
from config import Config
from compact import Compact
import string
import csv
from branch import Branch

def generate_instance():
    with open('inst_teste', 'w') as f:
        machine = random.randint(3, 5)
        jobs = random.randint(3, 5)
        f.write('{} {}\n'.format(jobs, machine))
        for j in range(jobs):
            order = list(range(machine))
            random.shuffle(order)
            time = random.randint(1, 10)
            str = '{} {}'.format(order[0], time)
            for i in range(1, machine):
                time = random.randint(1, 10)
                str += ' {} {}'.format(order[i], time)
            f.write('{}\n'.format(str))
        f.close()

def find_clique_applegate(compact: Compact):
    lo = compact.config.get_property('applegate_clique_lo')
    hi = compact.config.get_property('applegate_clique_hi')
    print(lo, hi)
    # input()
    constraints = 1
    total = 0
    while constraints > 0:
        compact.model.optimize(relax=True)
        constraints = 0
        hasCuts = 0
        for a in range(compact.instance.m):
            clique_cuts = compact.clique_cuts_best(a, lo, hi)
            constraints += clique_cuts
            total += clique_cuts
        print(constraints)
        # compact.model.write('teste.lp')
        compact.iterationsCuts += 1
        # input()
    return total

def find_general_clique(compact: Compact):
    selects = [0, 3]
    params = list(range(3, 9))
    for s in selects:
        for p in params:
            compact.config.conf['mip_general_cliques_select'] = s
            compact.config.conf['mip_general_cliques_parameter'] = p
            print(compact.config.conf)
            total = compact.mip_general_cliques()
            if total > 1:
                return total
    return total

def main(parameter):
    for i in range(100000):
        generate_instance()
        conf = Config('inst_teste.cfg')
        compact = Compact(conf)
        compact.instance.print()
        compact.constructProblemMcCormick()
        compact.model.verbose = 0
        compact.model.optimize()
        first_opt = compact.model.objective_value
        find_clique_applegate(compact)
        general_cliques = find_general_clique(compact)
        print('general', general_cliques)
        if general_cliques > 0:
            copyfile('inst_teste', 'teste{}'.format(i))
            compact.model.write('teste_clique_{}.lp'.format(i))
        compact.model.optimize()
        second_opt = compact.model.objective_value
        print(first_opt, second_opt)
        if abs(first_opt - second_opt) > 1e-8:
            copyfile('inst_teste', 'teste_erro{}'.format(i))
            compact.model.write('teste_clique_erro_{}.lp'.format(i))

        compact.printSolution()


    # print("teste")


if __name__ == "__main__":
    main()
