CXXFLAGS = -O3 -std=c++11 -mcmodel=medium
CXX = g++
LIBS = -lbifrost -pthread -lz

all: rnabridge-denovo clean

utils.o: utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

graph.o: graph.cpp 
	$(CXX) $(CXXFLAGS) -march=native -Wextra -c graph.cpp
bridging.o: bridging.cpp
	$(CXX) $(CXXFLAGS) -c bridging.cpp

rnabridge-denovo: graph.o bridging.o utils.o main.cpp
	$(CXX) $(CXXFLAGS) -march=native -Wextra -o rnabridge-denovo main.cpp graph.o bridging.o utils.o $(LIBS)
	
.PHONY:clean


clean:
	rm *.o
