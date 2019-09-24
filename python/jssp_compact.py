from sys import argv
import time
from JSSPInstance import JSSPInstance
from compact import Compact

inst = JSSPInstance(argv[1])

compact = Compact(inst)

# build the model
compact.constructProblem()

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
