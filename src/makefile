CXXFLAGS = -O3 -std=c++11 -mcmodel=medium
CC = g++

all: rnabridge-denova clean

utils.o: utils.cpp
	$(CC) $(CXXFLAGS) -c utils.cpp

graph.o: graph.cpp 
	$(CC) $(CXXFLAGS) -march=native -Wextra -c graph.cpp -lbifrost -pthread -lz

bridging.o: bridging.cpp
	$(CC) $(CXXFLAGS) -c bridging.cpp

rnabridge-denova: graph.o bridging.o utils.o main.cpp
	$(CC) $(CXXFLAGS) -march=native -Wextra -o rnabridge-denova main.cpp graph.o bridging.o utils.o -lbifrost -pthread -lz
	
.PHONY:clean


clean:
	rm *.o
