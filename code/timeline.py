time = input("Max Time")
jobs = input("Qtd jobs")
machines = input("Qtd machines")

timeline = [] 
vars = [] 

file = open("jssp.sol")

file.readline()

for line in file: 
  num,var,one,zero = line.split(" ")
  j,m0,t0,mf,tf = var.spplit(",") #x(1,i,1,3,1)
  x,j = j.split("(")
  vars.append((j,m0,t0,mf,tf))

