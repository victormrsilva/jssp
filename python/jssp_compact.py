from sys import argv
import time
from JSSPInstance import JSSPInstance
from compact import Compact
from branch import Branch

inst = JSSPInstance(argv[1])
# inst.print()
# input('instance')
compact = Compact(inst)

# build the model
# compact.constructProblemM()
# compact.constructProblemMcCormick()
# input('feito')
# compact.model.relax()
# compact.model.optimize()
# compact.printSolution()
# input()
# compact.constructProblemMcCormick()
compact.constructProblemMcCormickNonNegative()
# input()
branch = Branch()

branch.branch(compact)


start = time.time()
# execute with cutpool
# compact.optmizeCuts()

# # execute with adding cuts after relaxing the problem
compact.relax()
#
# # execute the integer problem
# compact.optimizeInteger()


end = time.time()

# printing results
# compact.printSolution()
print('Elapsed time: {}'.format(end - start))
