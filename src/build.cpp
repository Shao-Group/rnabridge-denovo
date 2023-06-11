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

void print_help()
{
	printf("\n");
	printf("Usage: rnabridge-denovo <input-read1> <input-read2> <output-bridge-sequence> <tmpfile-location> <threshold>\n");
	printf("\n");
}

graph* mygraph;
int main(int argc, char **argv)
{
	if(argc < 6)
	{
		printf("Parameters Error\n");
		print_help();
		return 0;
	}
	tmp = argv[4];
	mkdir(tmp.c_str(), 0755);

	mygraph = new graph(atoi(argv[5]));

	read(argv);
	mygraph->build(rd1, rd2, tmp);

	return 0;
}