CC=g++
OPTIONS=-Wall -g3 -DGRB -DDEBUG -std=c++11 -I/opt/gurobi752/linux64/include 
#-I/opt/gurobi752/linux64/include -L/opt/gurobi752/linux64/lib/ -lgurobi_c++ -lgurobi75 -lm `pkg-config --cflags cbc`
#-I/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/include/ `pkg-config --cflags cbc`

all:jssp

jssp:main.o Instance.o lp.o Compact.o 
	$(CC) $(OPTIONS) main.o Instance.o lp.o Compact.o -L/opt/gurobi752/linux64/lib -lgurobi75 -o jssp 
#-lasan -lm `pkg-config --libs cbc` 

main.o:main.cpp 
	$(CC) $(OPTIONS) -c main.cpp -o main.o

Instance.o:Instance.cpp Instance.hpp
	$(CC) $(OPTIONS) -c Instance.cpp -o Instance.o

lp.o:lp.cpp lp.hpp lp.h
	$(CC) $(OPTIONS) -c lp.cpp -o lp.o

Compact.o:Compact.cpp Compact.hpp Instance.o lp.o
	$(CC) $(OPTIONS) -c Compact.cpp -o Compact.o

clean:
	rm -rf *.o *~ jssp