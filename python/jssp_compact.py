from sys import argv
from JSSPInstance import JSSPInstance
from compact import Compact

inst = JSSPInstance(argv[1])

compact = Compact(inst)

# build the model
compact.constructProblem()

# execute with cutpool
# compact.optmizeCuts()

# # execute with adding cuts after relaxing the problem
compact.relax()
#
# # execute the integer problem
# compact.optimizeInteger()

# printing results
compact.printSolution()