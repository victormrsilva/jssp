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
compact.instance.print()
# input()
# print(conf.get_property("instance_name"))
# print(conf.get_property("instance"))
# input()
#
# if len(argv) < 3:
#     inst = JSSPInstance(argv[1], -1)
# else:
#     inst = JSSPInstance(argv[1], int(argv[2]))

instance_name = compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))

lo = compact.config.get_property('applegate_clique_lo')
hi = compact.config.get_property('applegate_clique_hi')
with open("{}.csv".format(instance_name), "a") as f:
    f.write("modo;iter;obj;time\n")

var = 'quad'
try:
    compact.iterationsCuts = 0
    compact.model.clear()

    compact.constructProblemM()
    compact.model.verbose = 0
    compact.model.optimize(relax=True)
    compact.config.conf['mip_general_cliques_quad'] = 1
    compact.config.conf['mip_general_cliques_not_sucessor'] = 0
    compact.config.conf['mip_general_cliques_sucessor'] = 0
    compact.config.conf['clique_cuts'] = 0
    compact.config.conf['mip_general_cliques_cuts'] = 0
    start = time.time()
    doitReturnValue = func_timeout(10800, compact.relax)
    # doitReturnValue = func_timeout(10800, compact.testCliqueMIP, args=(lo, h))
    end = time.time()
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
    print("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
except FunctionTimedOut:
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
    print("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
except Exception as e:
    print('error ', e)    

var = 'notsucessor'
try:
    compact.iterationsCuts = 0
    compact.model.clear()

    compact.constructProblemM()
    compact.model.verbose = 0
    compact.model.optimize(relax=True)
    basic_cuts = 0
    compact.config.conf['mip_general_cliques_quad'] = 0
    compact.config.conf['mip_general_cliques_not_sucessor'] = 1
    compact.config.conf['mip_general_cliques_sucessor'] = 0
    compact.config.conf['clique_cuts'] = 0
    compact.config.conf['mip_general_cliques_cuts'] = 0
    start = time.time()
    doitReturnValue = func_timeout(10800, compact.relax)
    # doitReturnValue = func_timeout(10800, compact.testCliqueMIP, args=(lo, h))
    end = time.time()
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
    print("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
except FunctionTimedOut:
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
    print("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
except Exception as e:
    print('error ', e)    

var = 'sucessor'
try:
    compact.iterationsCuts = 0
    compact.model.clear()

    compact.constructProblemM()
    compact.model.verbose = 0
    compact.model.optimize(relax=True)
    basic_cuts = 0
    compact.config.conf['mip_general_cliques_quad'] = 0
    compact.config.conf['mip_general_cliques_not_sucessor'] = 0
    compact.config.conf['mip_general_cliques_sucessor'] = 1
    compact.config.conf['clique_cuts'] = 0
    compact.config.conf['mip_general_cliques_cuts'] = 0
    start = time.time()
    doitReturnValue = func_timeout(10800, compact.relax)
    # doitReturnValue = func_timeout(10800, compact.testCliqueMIP, args=(lo, h))
    end = time.time()
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
    print("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
except FunctionTimedOut:
    compact.model.write('{}_{}.lp'.format(instance_name, var))
    compact.model.optimize(relax=True)
    with open("{}.csv".format(instance_name), "a") as f:
        f.write("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
    print("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
except Exception as e:
    print('error ', e)    


# basic + general
for i in range(5, 8):
    var = 'basic+general{}'.format(i)
    try:
        compact.iterationsCuts = 0
        compact.model.clear()

        compact.constructProblemM()

        compact.model.verbose = 0
        compact.model.optimize(relax=True)
        basic_cuts = 0
        compact.config.conf['clique_cuts'] = 0
        compact.config.conf['mip_general_cliques_sucessor'] = 1
        compact.config.conf['mip_general_cliques_not_sucessor'] = 1
        compact.config.conf['mip_general_cliques_quad'] = 1
        compact.config.conf['mip_general_cliques_cuts'] = 1
        compact.config.conf['mip_general_cliques_parameter'] = i
        start = time.time()
        doitReturnValue = func_timeout(10800, compact.relax)
        # doitReturnValue = func_timeout(10800, compact.testCliqueMIP, args=(lo, h))
        end = time.time()
        compact.model.write('{}_{}.lp'.format(instance_name, var))
        compact.model.optimize(relax=True)
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
        print("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
    except FunctionTimedOut:
        compact.model.write('{}_{}.lp'.format(instance_name, var))
        compact.model.optimize(relax=True)
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
        print("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
    except Exception as e:
        print('error ', e)

# basic + applegate + general
for i in range(5, 8):
    var = 'basic+applegate+general{}'.format(i)
    try:
        compact.iterationsCuts = 0
        compact.model.clear()

        compact.constructProblemM()

        compact.model.verbose = 0
        compact.model.optimize(relax=True)
        
        compact.config.conf['clique_cuts'] = 1
        compact.config.conf['mip_general_cliques_sucessor'] = 1
        compact.config.conf['mip_general_cliques_not_sucessor'] = 1
        compact.config.conf['mip_general_cliques_quad'] = 1
        compact.config.conf['mip_general_cliques_cuts'] = 1
        compact.config.conf['mip_general_cliques_parameter'] = i
        start = time.time()
        doitReturnValue = func_timeout(10800, compact.relax)
        # doitReturnValue = func_timeout(10800, compact.testCliqueMIP, args=(lo, h))
        end = time.time()
        compact.model.write('{}_{}.lp'.format(instance_name, var))
        compact.model.optimize(relax=True)
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
        print("{};{};{};{}\n".format(var, compact.iterationsCuts, compact.model.objective_value, end - start))
    except FunctionTimedOut:
        compact.model.write('{}_{}.lp'.format(instance_name, var))
        compact.model.optimize(relax=True)
        with open("{}.csv".format(instance_name), "a") as f:
            f.write("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
        print("{};{};{};timeout\n".format(var, compact.iterationsCuts, compact.model.objective_value))
    except Exception as e:
        print('error ', e)


