CC=g++
OPTIONS=-Wall -g3  -DGRB -std=c++11 -I/opt/gurobi752/linux64/include 
#-I/opt/gurobi752/linux64/include -L/opt/gurobi752/linux64/lib/ -lgurobi_c++ -lgurobi75 -lm `pkg-config --cflags cbc`
#-I/usr/lib/gcc/x86_64-pc-linux-gnu/7.2.0/include/ `pkg-config --cflags cbc`

all:jssp

jssp:main.o Instance.o lp.o Flow.o Compact.o Kondili.o
	$(CC) $(OPTIONS) main.o Instance.o lp.o Flow.o Compact.o Kondili.o -L/opt/gurobi800/linux64/lib -lgurobi80 -o jssp 
#-lasan -lm `pkg-config --libs cbc` 

main.o:main.cpp 
	$(CC) $(OPTIONS) -c main.cpp -o main.o

Instance.o:src/Instance.cpp src/Instance.hpp
	$(CC) $(OPTIONS) -c src/Instance.cpp -o Instance.o

lp.o:src/lp.cpp src/lp.hpp src/lp.h
	$(CC) $(OPTIONS) -c src/lp.cpp -o lp.o

Compact.o:src/Compact.cpp src/Compact.hpp Instance.o lp.o
	$(CC) $(OPTIONS) -c src/Compact.cpp -o Compact.o

Kondili.o:src/Kondili.cpp src/Kondili.hpp Instance.o lp.o
	$(CC) $(OPTIONS) -c src/Kondili.cpp -o Kondili.o

Flow.o:src/Flow.cpp src/Flow.hpp Instance.o lp.o
	$(CC) $(OPTIONS) -c src/Flow.cpp -o Flow.o

clean:
	rm -rf *.o *~ jssp