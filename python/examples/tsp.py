#!/usr/bin/python
# Utilizacao: gurobi.sh tsp.py instancia.txt

import sys
import math
import itertools
from gurobipy import *

n = 0
EPS = 0.0001

def main():
    if len(sys.argv) < 1:
        print('\n  Indique o arquivo de entrada. Uso correto: gurobi.sh tsp.py instancia.txt\n')
        return

    # lendo coordenadas do arquivo de entrada
    problem_file = open(sys.argv[1], 'r')
    n = int(problem_file.readline())
    points = []
    for i in range(n):
        line = problem_file.readline().split()
        points.append((int(line[0]), int(line[1])))

    # calculando distancias euclidianas
    dist = {(i,j) :
        math.sqrt(sum((points[i][k]-points[j][k])**2 for k in range(2)))
        for i in range(n) for j in range(n)}
    # criando o modelo inicial
    model = create_model(n, points, dist)
    while True:
        model.optimize()
        #print model.getVars()
        print('Custo da solucao relaxada: %.5f' % model.objVal)
       
        
        if not add_cut(model):
            break

def create_model(n, points, dist):
    model = Model()
    # crie o modelo aqui...
    x = [[ 0 for x in range(n)] for y in range(n)]
    for i in range(n):
      for j in range(n):
        if (i == j):
          continue
        x[i][j] = model.addVar(lb=0,ub=1,obj=dist[i,j], name="x_(%d,%d)" % (i,j))
        
    for i in range(n):
      model.addConstr(sum([ x[i][j] for j in range(n) if i != j]) == 1, name="x_%d,j" % (i))

    for j in range(n):
      model.addConstr(sum([ x[i][j] for i in range(n) if i != j ]) == 1, name="x_i,%d" % (j))
    
    model._n = n
    return model

def add_cut(model):
    global EPS
    edges = []
    for i in range(model._n):
      for j in range(model._n):
        if (i == j):
          continue
        if model.getVarByName("x_(%d,%d)" % (i,j)).X > EPS:
          edges.append((i,j))
    #print edges
    
    tour = subtour(edges,model._n)
    #print tour
    
    if len(tour) < model._n:
      model.addConstr(sum(model.getVarByName("x_(%d,%d)" % (i,j)) for (i,j) in tour) <= len(tour) - 1 + EPS);
      return True
      
        
    return False

def subtour(edges,n):
    """
    Dada uma lista de arestas conectadas, encontra a menor sub-rota
    (extraido dos exemplos do gurobi)
    """
    edges = tuplelist(edges)
    unvisited = list(range(n))
    cycle = range(n+1) # initial length has 1 more city
    print(edges)
    input()

    while unvisited: # true if list is non-empty
        thiscycle = []
        neighbors = unvisited
        print('neigh', neighbors)
        while neighbors:
            current = neighbors[0]
            print('cur', current)
            thiscycle.append(current)
            print('this', thiscycle)
            unvisited.remove(current)
            print('unvisited', unvisited)
            neighbors = [j for i,j in edges.select(current,'*') if j in unvisited]
            print('neigh', neighbors)
            input()
        if len(cycle) > len(thiscycle) and len(thiscycle) > 1:
            print('cycle', cycle, 'this', thiscycle)
            cycle = thiscycle

    result = []
    prev = cycle[0]
    for i in range(1,len(cycle)):
        result.append((prev,cycle[i]))
        prev = cycle[i]
    result.append((prev,cycle[0]))
    return result

if __name__ == "__main__":
    main()