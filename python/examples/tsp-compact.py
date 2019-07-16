from sys import argv
from typing import List, Tuple
import networkx as nx
from tspdata import TSPData
from mip.model import Model, xsum, BINARY
from mip.callbacks import CutsGenerator, CutPool
from itertools import product


class SubTourCutGenerator(CutsGenerator):
    def __init__(self, Fl: List[Tuple[int, int]], m: Model):
        self.F = Fl
        self.m = m

    def generate_cuts(self, model: Model):
        G = nx.DiGraph()
        r = [(v, v.x) for v in model.vars if v.name.startswith('x(')]
        U = [int(v.name.split('(')[1].split(',')[0]) for v, f in r]
        V = [int(v.name.split(')')[0].split(',')[1]) for v, f in r]

        for v in model.vars:
            print('{} = {} '.format(v.name, v.x), end='')
        print()
        input()
        cp = CutPool()
        for i in range(len(U)):
            G.add_edge(U[i], V[i], capacity=r[i][1])
        for (u, v) in F:
            if u not in U or v not in V:
                continue
            val, (S, NS) = nx.minimum_cut(G, u, v)
            if val <= 0.99:
                arcsInS = [(v, f) for i, (v, f) in enumerate(r)
                           if U[i] in S and V[i] in S]
                if sum(f for v, f in arcsInS) >= (len(S)-1)+1e-4:
                    cut = xsum(1.0*v for v, fm in arcsInS) <= len(S)-1
                    cp.add(cut)
                    if len(cp.cuts) > 256:
                        for cut in cp.cuts:
                            model.add_cut(cut)
                        return
        for cut in cp.cuts:
            model.add_cut(cut)
        return


inst = TSPData(argv[1])
n, d = inst.n, inst.d

model = Model()

x = [[model.add_var(name='x({},{})'.format(i, j),
                    var_type=BINARY) for j in range(n)] for i in range(n)]
y = [model.add_var(name='y({})'.format(i),
                   lb=0.0, ub=n) for i in range(n)]

model.objective = xsum(d[i][j] * x[i][j] for j in range(n) for i in range(n))

for i in range(n):
    model += xsum(x[j][i] for j in range(n) if j != i) == 1
    model += xsum(x[i][j] for j in range(n) if j != i) == 1
for (i, j) in [(i, j) for (i, j) in
               product(range(1, n), range(1, n)) if i != j]:
    model += y[i] - (n + 1) * x[i][j] >= y[j] - n

F = []
for i in range(n):
    (md, dp) = (0, -1)
    for j in [k for k in range(n) if k != i]:
        if d[i][j] > md:
            (md, dp) = (d[i][j], j)
    F.append((i, dp))

model.cuts_generator = SubTourCutGenerator(F, model)
model.optimize()

arcs = [(i, j) for i in range(n) for j in range(n) if x[i][j].x >= 0.99]
print('optimal route : {}'.format(arcs))