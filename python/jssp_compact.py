from sys import argv
import time
from JSSPInstance import JSSPInstance
from compact import Compact
import string
import csv
from branch import Branch

inst = JSSPInstance(argv[1])
# inst.print()
# input('instance')
compact = Compact(inst)
instance_name = compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))

# build the model
# compact.constructProblemM()
# compact.constructProblemMcCormick()
# input('feito')
# compact.model.relax()
# compact.model.optimize()
# compact.printSolution()
# input()
# compact.constructProblemMcCormick()

# compact.constructProblemMcCormickNonNegative()
# input()
# branch = Branch()

# branch.branch(compact)


# start = time.time()
# execute with cutpool
# compact.optmizeCuts()

# # execute with adding cuts after relaxing the problem
file = open('{}.csv'.format(instance_name), "w")
writer = csv.writer(file)
# writer.writerow(['type', 'lc', 'hc', 'first', 'last', 'time', 'cuts'])
# maxlc = min(3, compact.instance.n)
# maxhc = min(8, compact.instance.n+1)
# for lc in range(2, maxlc):
#     for hc in range(3, maxhc):
#         print('Executing {} MIP with lc = {} and hc = {} ...'.format(instance_name, lc, hc))
#         compact = Compact(inst)
#         compact.constructProblemMcCormickNonNegative()
#         first, last, time, cuts = compact.testCliqueMIP(lc, hc)
#         writer.writerow(['mip', lc, hc, first, last, time, cuts])
#         compact.model.write('{}_clique_MIP_lc_{}_hc_{}.lp'.format(instance_name, lc, hc))
#         print(' Done!!!')
#         del compact

writer.writerow(['type', 'try', 'maxstep', 'first', 'last', 'time', 'cuts'])
maxsteps = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
for step in maxsteps:
    for i in range(5):
        del compact
        print('Executing {} SA with maxsteps = {} and iteration = {} ...'.format(instance_name, step, i))
        compact = Compact(inst)
        compact.constructProblemMcCormickNonNegative()
        first, last, time, cuts = compact.testCliqueSA(step)
        writer.writerow(['SA', i, step, first, last, time, cuts])
        compact.model.write('{}_clique_SA_{}_iter_{}.lp'.format(instance_name, step, i))
        print(' Done!!!')
file.close()

# print(compact.testCliqueMIP(2, 5))
# print(compact.testCliqueSA())
# compact.model.write(
#     '{}_modelMCLinNonNeg_cuts.lp'.format(compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))))

#
# # execute the integer problem
# compact.optimizeInteger()


# end = time.time()

# printing results
# compact.printSolution()
# print('Elapsed time: {}'.format(end - start))
