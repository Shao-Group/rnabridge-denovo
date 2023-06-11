CXXFLAGS = -O3 -std=c++11 -mcmodel=medium -g
CXX = g++
LIBS = -lbifrost -pthread -lz

all: rnabridge-denovo build clean

utils.o: utils.cpp
	$(CXX) $(CXXFLAGS) -c utils.cpp

graph.o: graph.cpp 
	$(CXX) $(CXXFLAGS) -march=native -Wextra -c graph.cpp
bridging.o: bridging.cpp
	$(CXX) $(CXXFLAGS) -c bridging.cpp

build: graph.o bridging.o utils.o main.cpp
	$(CXX) $(CXXFLAGS) -march=native -Wextra -o rnabridge-build build.cpp graph.o bridging.o utils.o $(LIBS) ${LDFLAGS}

rnabridge-denovo: graph.o bridging.o utils.o main.cpp
	$(CXX) $(CXXFLAGS) -march=native -Wextra -o rnabridge-denovo main.cpp graph.o bridging.o utils.o $(LIBS) ${LDFLAGS}
	
.PHONY:clean


clean:
	rm *.o
