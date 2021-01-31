CXXFLAGS = -O3 -std=c++11 -mcmodel=medium
CC = g++

all: find_bridge graph_align bridge2fa

find_bridge: find_bridge.cpp
	$(CC) $(CXXFLAGS) -o find_bridge find_bridge.cpp

bridge2fa: bridge2fa.cpp
	$(CC) $(CXXFLAGS) -o bridge2fa bridge2fa.cpp

graph_align: graph_align.cpp
	$(CC) $(CXXFLAGS) -march=native -lbifrost -pthread -lz -Wextra -o graph_align graph_align.cpp