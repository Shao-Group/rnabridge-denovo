#include <iostream>
#include <stack>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "bifrost/CompactedDBG.hpp"
#include "bifrost/ColoredCDBG.hpp"
#include "graph.h"
#include "bridging.h"

using namespace std;

string rd1, rd2, output, tmp;

void read(char **argv)
{
	rd1 = argv[1];
	rd2 = argv[2];
	output = argv[3];
}

graph mygraph;
bridging findbridge;

int main(int argc, char **argv)
{
	if(argc < 4)
	{
		printf("Parameters Error\n");
		return 0;
	}
	tmp = "./tmpfiles/";
	mkdir(tmp.c_str(), 0755);

	read(argv);
	mygraph.build(rd1, rd2, tmp);
	findbridge.read(tmp + "/graph", tmp + "query");
	findbridge.findPaths();
	findbridge.write(output);

	rmdir(tmp.c_str());
	return 0;
}