from sys import argv
import time
from JSSPInstance import JSSPInstance
from config import Config
from compact import Compact
import string
import csv
from branch import Branch
from func_timeout import func_timeout, FunctionTimedOut


start = time.time()
conf = Config(argv[1])
compact = Compact(conf)
# compact.instance.print()
# input()
# print(conf.get_property("instance_name"))
# print(conf.get_property("instance"))
# input()
#
# if len(argv) < 3:
#     inst = JSSPInstance(argv[1], -1)
# else:
#     inst = JSSPInstance(argv[1], int(argv[2]))

# inst.print()
# input('instance')
# compact = Compact(inst)
instance_name = compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))

# build the model
# compact.constructProblemMSubCycle()
# compact.model.relax()
# compact.model.optimize()
# compact.triangle_cuts_best(0)
# input('teste')
# compact.optimizeSubCycle()
# compact.model.optimize()
# compact.printSolution()
# print(compact.instance.o)
# input()
# compact.constructProblemMcCormickNonNegative()

# compact.noname_clique()

# str = argv[1]

# str = str + ";{0:.6f}".format(compact.model.objective_value)


lo = compact.config.get_property('applegate_clique_lo')
hi = compact.config.get_property('applegate_clique_hi')
hi =  min(hi, compact.instance.n)
file = open('{}.csv'.format(instance_name), "w")
with open("{}.csv".format(instance_name), "a") as f:
    f.write("lo;hi;obj;time\n")
for i in range(hi, 5, -1):
    try:
        if conf.get_property("problem") == 0:  # Big-M
            compact.constructProblemM()
        elif conf.get_property("problem") == 1:  # McCormick
            compact.constructProblemMcCormick()
        elif conf.get_property("problem") == 2:  # McCormick Non Negative
            compact.constructProblemMcCormickNonNegative()
        compact.model.verbose = 0
        compact.model.optimize(relax=True)
        start = time.time()
        doitReturnValue = func_timeout(10, compact.testCliqueMIP, args=(lo, i))
        end = time.time()
        print(doitReturnValue)
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};{}\n".format(lo, i, compact.model.objective_value, end - start))
        print(lo, i, compact.model.objective_value, end - start)
    except FunctionTimedOut:
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};-\n".format(lo, i, compact.model.objective_value))
        print(lo, i, compact.model.objective_value, 'timeout')
    # except Exception as e:
        # print('test')
# compact.testCliqueMIP(lo, hi)
# compact.printSolution()
# input()
# compact.mip_general_cliques()
# compact.printSolution()
# input()
# compact.model.optimize(relax=True)
# str = str + ";{0:.6f}".format(compact.model.objective_value)
# compact.model.optimize()
# compact.printSolution()
# input()
# str = str + ";{0:.6f}".format(compact.model.objective_value)
# end = time.time()
# str = str + ";{0:.6f}".format(end-start)

# with open("{}.csv".format(instance_name), "a") as f:
#     f.write(str+"\n")
# compact.testCliqueMIP(2, 4)
# input('feito')
# compact.model.relax()
# compact.model.optimize()
# compact.printSolution()
# compact.model.write('teste.lp')
# input('teste')
# # compact.constructProblemMcCormick()
#
# # compact.constructProblemMcCormickNonNegative()
# # input()
# # branch = Branch()
#
# # branch.branch(compact)
#
#
# # start = time.time()
# # execute with cutpool
# # compact.optmizeCuts()
#
# # # execute with adding cuts after relaxing the problem
# file = open('{}.csv'.format(instance_name), "w")
# writer = csv.writer(file)
# # writer.writerow(['type', 'lc', 'hc', 'first', 'last', 'time', 'cuts'])
# # maxlc = min(3, compact.instance.n)
# # maxhc = min(8, compact.instance.n+1)
# # for lc in range(2, maxlc):
# #     for hc in range(3, maxhc):
# #         print('Executing {} MIP with lc = {} and hc = {} ...'.format(instance_name, lc, hc))
# #         compact = Compact(inst)
# #         compact.constructProblemMcCormickNonNegative()
# #         first, last, time, cuts = compact.testCliqueMIP(lc, hc)
# #         writer.writerow(['mip', lc, hc, first, last, time, cuts])
# #         compact.model.write('{}_clique_MIP_lc_{}_hc_{}.lp'.format(instance_name, lc, hc))
# #         print(' Done!!!')
# #         del compact
#
# writer.writerow(['type', 'try', 'maxstep', 'minimum', 'maximum', 'exact', 'first', 'last', 'time', 'cuts'])
# maxsteps = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
# for step in maxsteps:
#     for i in range(5):
#         del compact
#         print('Executing {} SA with maxsteps = {} and iteration = {}, with prob = [0.25 0.25 0.25 0.25] ...'.format(instance_name, step, i))
#         compact = Compact(inst)
#         compact.constructProblemMcCormickNonNegative()
#         first, last, time, cuts, minimum, maximum, exact = compact.testCliqueSA(step)
#         writer.writerow(['SA', i, step, minimum, maximum, exact, first, last, time, cuts])
#         # compact.model.write('{}_clique_SA_{}_iter_{}.lp'.format(instance_name, step, i))
#         print('SA iteration: {}, maxstep: {}, minimum: {}, maximum: {}, exact: {}, first_lp: {}'
#               ', last_lp: {}, time: {}, cuts: {}'.format(i, step, minimum, maximum, exact, first, last, time, cuts))
#         input()
# # writer.writerow(['type', 'try', 'maxstep', 'l', 'first', 'last', 'time', 'cuts'])
# # maxsteps = [1000, 5000, 10000, 50000, 100000, 500000, 1000000]
# # ls = [10, 25, 50]
# # for step in maxsteps:
# #     for l in ls:
# #         for i in range(5):
# #             del compact
# #             print('Executing {} LAHC with maxsteps = {} and iteration = {} ...'.format(instance_name, step, i))
# #             compact = Compact(inst)
# #             compact.constructProblemMcCormickNonNegative()
# #             first, last, time, cuts = compact.testCliqueLAHC(step, l)
# #             writer.writerow(['LAHC', i, step, l, first, last, time, cuts])
# #             compact.model.write('{}_clique_LAHC_{}_l_{}_iter_{}.lp'.format(instance_name, step, l, i))
#
# file.close()
#
# # print(compact.testCliqueMIP(2, 5))
# # print(compact.testCliqueSA())
# # compact.model.write(
# #     '{}_modelMCLinNonNeg_cuts.lp'.format(compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))))
#
# #
# # # execute the integer problem
# # compact.optimizeInteger()
#
#
# # end = time.time()
#
# # printing results
# # compact.printSolution()
# # print('Elapsed time: {}'.format(end - start))