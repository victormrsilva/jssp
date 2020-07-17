from JSSPInstance import JSSPInstance
from sys import argv
import string


import numpy as np




def solve(b, x, t):
    print(x, t, b)
    m = Model(solver_name="grb")
    m.verbose = 0
    a0 = m.add_var(var_type=INTEGER, lb=1, name="b")
    a = [m.add_var(var_type=INTEGER, lb=1, name="a{}".format(i)) for i in range(len(x))]
    print(a0, a)
    m.objective = a0
    for i in range(len(x)):
        m += a[i] - t[i]*a0 <= 0, 'restr{}'.format(i)
    m += xsum(x[i]*a[i] for i in range(len(x))) - b*a0 >= 0, 'obj'

    m.write('teste.lp')
    m.optimize()
    return a0.x, [a[i].x for i in range(len(x))]

def teste():
    v = np.array([[0, 4, 15], [0, 18, 15], [0, 4, 15], [18, 4, 15], [18, 18, 15], [18, 18, 15]])
    for i in range(len(v)):
        for k in range(i+1, len(v)):
            for l in range(k+1, len(v)):
                aux = v[[i,k,l],:]
                a0 = round(np.linalg.det(aux), 1)
                str = ''
                for j in range(3):
                    aux2 = aux.copy()
                    aux2[:,j] = 1
                    a = round(np.linalg.det(aux2), 1)
                    str += ' + {}*t{}'.format(a, j)
                str += ' >= {}'.format(a0)
                print(str)
    input()



if __name__ == "__main__":
    teste()
    conf = Config(argv[1])
    compact = Compact(conf)
    instance_name = compact.instance.instancename.translate(str.maketrans('', '', string.punctuation))

    out = open("{}_restr_res.csv".format(instance_name), "w")

    with open("{}_restr.csv".format(instance_name), "r") as f:
        line = f.readline()
        cnt = 1
        while line:
            print(line)
            line = f.readline()
            aux = line.split(';')
            size = len(aux)
            restr = aux[0]
            b = float(aux[2])
            name = list()
            x = list()
            t = list()

            i = 3
            if b > 0:
                while i < len(aux):
                    print(i, len(aux))
                    name.append(aux[i])
                    x.append(float(aux[i+1]))
                    t.append(float(aux[i+3]))
                    i = i + 4

                int_b, int_x = solve(b, x, t)
                
                print(int_b, int_x)
                tup = list()
                print(name)
                for i in range(len(name)):
                    test = name[i].split('(')[1].split(')')[0].split(',')
                    tup.append((int(test[0].strip('j')), int(test[1].strip('m'))))
                text = '{};{};{}'.format(restr, 'b', int_b)
                for i in range(len(name)):
                    text += ';{};{};p;{};est;{}'.format(name[i], int_x[i], compact.instance.times[tup[i][0]][tup[i][1]], compact.instance.est[tup[i][0]][tup[i][1]])
                text += '\n'
                out.write(text)
                # input(text)