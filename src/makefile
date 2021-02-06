CXXFLAGS = -O3 -std=c++11 -mcmodel=medium
CC = g++
LDFLAGS = -lbifrost -pthread -lz

all: rnabridge-denovo clean

utils.o: utils.cpp
	$(CC) $(CXXFLAGS) -c utils.cpp

graph.o: graph.cpp 
	$(CC) $(CXXFLAGS) -march=native -Wextra -c graph.cpp
bridging.o: bridging.cpp
	$(CC) $(CXXFLAGS) -c bridging.cpp

rnabridge-denovo: graph.o bridging.o utils.o main.cpp
	$(CC) $(CXXFLAGS) -march=native -Wextra -o rnabridge-denovo main.cpp graph.o bridging.o utils.o $(LDFLAGS)
	
.PHONY:clean


clean:
	rm *.o
