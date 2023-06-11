#include <iostream>
#include <stack>
#include <algorithm>
#include <queue>
#include <map>
#include <set>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <unordered_set>
#include <sstream>
#include "utils.h"

using namespace std;

#ifndef __BRIDGING_H__
#define __BRIDGING_H__
const int kbottle = 1;

class nodeDistance{

public:
    int id;
    int pre;
    int dis[kbottle];
    int seq_len;
    int length;
    

    nodeDistance(int i, int p, int d[], int l, int sl)
    {
        id = i;
        pre = p;
        length = l;
        seq_len = sl;
        for(int i = 0; i < kbottle; i++)
            dis[i] = d[i];

    }

    bool operator < (const struct nodeDistance &n1) const
    {
        for(int i = 0; i < kbottle; i++)
            if(dis[i] != n1.dis[i])
                return dis[i] < n1.dis[i];
        return 1;
    }
};

struct readinfo{int s, t, id; bool found = 0; string l,r,path;};

class bridging{

    const int INF = 0x7fffffff/2;
    const static int maxN = 30000000;
    const static int maxn = 30000000;
    const int maxlen = 400;

    int n, m;
    vector<string> unis;
    vector<int> value;
    vector<int> unilen;

    int h[maxn<<1] = {0}, edgecnt = 0;
    struct edgelist{int to, next;} q[maxn<<2];

    int vst[maxn<<1], prevex[maxn<<1];
    int dst[maxn<<1][kbottle];
    int inpq[maxn<<1];

    struct dfs_stuck{int u,ch,length,seq_len;}sp[maxn<<1];
    int dvst[maxn<<1] = {0};

    vector<int> leaves;

    vector<readinfo> all[maxn<<1];
    int in[maxn<<1] = {0}, finish[maxn<<1] = {0};
    queue<int> s_queue;

    int undo = 0;
    vector<readinfo> ans;

    void add(int x, int y);
    string path2String(int s, int t, int k, string l, string r);
    vector<readinfo> dijkstra(int start, vector<readinfo> v);
    vector<readinfo> partlyDijkstra(int start, int pres, vector<readinfo> v);
    void release(int x);
    int findNext(int s);

public:
    void read(string graphinput, string queryinput);
    void findPaths();
    void write(string path);

};
#endif