#include <iostream>
#include <stack>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <fstream>
#include <unordered_set>
#include <sstream>

using namespace std;

string ans;

int main(int argc, char *argv[])
{
    ans = argv[1];
    FILE * readout = fopen (argv[2], "w");

    int cnt = 0;

    string line;

    ifstream fin1(ans);
    while(getline(fin1, line))
    {
        stringstream liness(line);

        string id;
        string bridge;

        getline(liness, id, '\t');
        getline(liness, bridge, '\t');
        cnt++;
        id = '>' + id;
        fprintf(readout, "%s\n%s\n", id.c_str(), bridge.c_str());
    }

    printf("Bridge: %d\n", cnt);
    return 0;
}
