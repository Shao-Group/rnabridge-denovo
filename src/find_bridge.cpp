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
#include <ctime>
#include <time.h>

using namespace std;
const int INF = 0x7fffffff/2;
const int maxN = 20000000;
const int maxn = 6000000;
const int maxlen = 400;
const int kbottle = 5;

int n, m;
vector<string> unis;
vector<int> value;
vector<int> unilen;

int h[maxn<<1] = {0}, edgecnt = 0;
struct edgelist{int to, next;}
q[maxn<<2];
int vst[maxn<<1], prevex[maxn<<1];
int dst[maxn<<1][kbottle];
int inpq[maxn<<1];

void Add(int x, int y)
{
    q[++edgecnt] = (edgelist){y, h[x]};
    h[x] = edgecnt;
}

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
int vnodes = 0;
vector<int> leaves;

string trans(string a)
{
    string b = "";
    for(int i = a.length() - 1; i >= 0; i--)
    {
        if(a[i] == 'A')
            b = b + 'T';
        if(a[i] == 'G')
            b = b + 'C';
        if(a[i] == 'C')
            b = b + 'G';
        if(a[i] == 'T')
            b = b + 'A';
    }

    return b;
}

string path2string(int s, int t, int k, string l, string r)
{
    int x = t;
    string tmp;

    int len = l.length();
    string ret = "";

    string lt = l.substr(len - k);
    string rt = r.substr(0, k);

    if(s == t)
    {
        tmp = unis[x>>1];
        if(x&1)
            tmp = trans(tmp);

        size_t ls = tmp.find(l);
        size_t rs = tmp.find(r);            
        if(ls != string::npos && rs != string::npos)
        {
            rs = rs + len;
            return tmp.substr(ls, rs - ls);
        }

        if(ls != string::npos)
        {
            rs = tmp.find(rt);
            return tmp.substr(ls, rs - ls) + r;
        }
        
        if(rs != string::npos)
        {
            rs = rs + len;
            ls = tmp.find(lt) + k;
            return l + tmp.substr(ls, rs - ls);
        }

        ls = tmp.find(lt) + k;
        rs = tmp.find(rt);

        if(rs < ls)
        {
            if(ls - rs > r.length())
                return "-1";
            return l + r.substr(ls - rs);
        }
        return l + tmp.substr(ls, rs - ls) + r;
    }


    // if(prevex[t] == s)
    // {
    //     int y = prevex[x];
    //     tmp = unis[x>>1];
    //     string tmpy = unis[y>>1];
    //     if(x&1)
    //         tmp = trans(tmp);
    //     if(y&1)
    //         tmpy = trans(tmpy);
    //     int ls = tmpy.find(lt) + k;
    //     int rs = tmp.find(rt);

    //     string ll = l + tmpy.substr(ls);
    //     string rr = tmp.substr(0, rs) + r;

    //     return ll + rr.substr(k-1);

    // }


    tmp = unis[x>>1];
    if(x&1)
        tmp = trans(tmp);
    int rs = tmp.find(rt);
    ret = tmp.substr(0, rs) + r;
    x = prevex[x];

    while(x != s)
    {
        tmp = unis[x>>1];
        if(x&1)
            tmp = trans(tmp);
        // if(s == 421051 && t == 420668)
        //     cout<<x<<" "<<prevex[x]<<endl;
        //if(s == 277394 && t== 277322)
            //cout<<x<<" "<<prevex[x]<<endl;
        ret = tmp.substr(0, tmp.length() - k + 1) + ret;
        //ret = tmp + ret;
        x = prevex[x];
    }

    tmp = unis[x>>1];
    if(x&1)
        tmp = trans(tmp);
    rs = tmp.find(lt) + k;

    ret = l + tmp.substr(rs) + ret.substr(k-1);

    return ret;
}

vector<readinfo> Djikstra(int start, vector<readinfo> v)
{
    vector<readinfo> ret;
    if(v.size() == 0)
        return ret;

    int sdis[kbottle];
    for(int i = 0; i < kbottle; i++)
        sdis[i] = INF;

    nodeDistance source(start, -1, sdis, 0, unilen[start>>1]);
    for(int i = 0; i < n<<1 ; i++)
    {
        for(int j = 0; j < kbottle; j++)
            dst[i][j] = -INF;
        prevex[i] = -1;
    }

    leaves.clear();

    priority_queue<nodeDistance> pq;
    pq.push(source);

    nodeDistance n(-1, -1, sdis, 0, 0);

    int x, num = 0;

    while(pq.size() > 0)
    {
        n = pq.top();
        pq.pop();
        
        x = n.id;

        while(vst[x] == start && !pq.empty())
        {
            n = pq.top();
            pq.pop();
            
            x = n.id;            
        }

        if(vst[x] == start && pq.empty())
            break;


        vst[x] = start;
        prevex[x] = n.pre;
        memcpy(dst[x], n.dis, sizeof(n.dis));
        num++;
        vnodes++;


        if(n.seq_len > maxlen)
        {
            leaves.push_back(x);
            continue;
        }

        for (int i = h[x];i;i = q[i].next)
        {   
            int y = q[i].to;
            int tmpd[kbottle];

            for(int j = 0; j < kbottle; j++)
            {
                if(value[y] < dst[x][j])
                {
                    tmpd[j] = value[y];
                    for(int k = j + 1; k < kbottle; k++)
                        tmpd[k] = dst[x][k-1];
                    break;
                }

                tmpd[j] = dst[x][j];
            }

            if(vst[y] != start)
                for(int j = 0; j < kbottle; j++)
                    if(dst[y][j] < tmpd[j])
                    {
                        memcpy(dst[y], tmpd, sizeof(tmpd));
                        pq.push(nodeDistance(y, x, dst[y], n.length + 1, n.seq_len + unilen[y>>1] - 30));
                        break;
                    }
                    else if(dst[y][j] > tmpd[j])
                        break;
        }
    }

    for(int i = 0; i < v.size(); i++)
        if(vst[v[i].t] == start)
        {
            v[i].found = true;
            v[i].path = path2string(v[i].s, v[i].t, 31, v[i].l, v[i].r);

            if(v[i].path != "-1")
                ret.push_back(v[i]);
        }

    return ret;
}

struct dfs_stuck{int u,ch,length,seq_len;}sp[10000001];

int dvst[maxn<<1] = {0};

vector<readinfo> PartlyDjikstra(int start, int pres, vector<readinfo> v)
{
    int tp = 1;
    priority_queue<nodeDistance> pq;

    for(int i = 0; i < kbottle; i++)
        dst[start][i] = INF;

    prevex[start] = -1;

    sp[tp].u = start;
    sp[tp].ch = h[start];
    sp[tp].length = 0;
    sp[tp].seq_len = 0;
    int u, w;

    while(tp)
    {
        dfs_stuck &cur = sp[tp];
        u = cur.u;
        vst[u] = start;

        //if(cur.ch == h[u] && start == 78)
            //printf("UP %d %d %d %d %d\n", u, dst[u], prevex[u], value[u], cur.length);
        if(cur.ch == 0)
        {
            tp--;
            continue;
        }

        w = q[cur.ch].to;
        cur.ch = q[cur.ch].next;

        if(vst[w] == pres && prevex[w] == u)
        {           

            for(int j = 0; j < kbottle; j++)
            {
                if(value[w] < dst[u][j])
                {
                    dst[w][j] = value[w];
                    for(int k = j + 1; k < kbottle; k++)
                        dst[w][k] = dst[u][k-1];
                    break;
                }

                dst[w][j] = dst[u][j];
            }

            tp++;
            sp[tp].u = w;
            sp[tp].ch = h[w];
            sp[tp].length = cur.length + 1;
            sp[tp].seq_len = cur.seq_len + unilen[w>>1] - 30;
        }
        else
        {            
            int tmpd[kbottle];

            for(int j = 0; j < kbottle; j++)
            {
                if(value[w] < dst[u][j])
                {
                    tmpd[j] = value[w];
                    for(int k = j + 1; k < kbottle; k++)
                        tmpd[k] = dst[u][k-1];
                    break;
                }

                tmpd[j] = dst[u][j];
            }


            if(inpq[w] != start)
            {
                memcpy(dst[w], tmpd, sizeof(tmpd));
                pq.push(nodeDistance(w, u, dst[w], cur.length + 1, cur.seq_len + unilen[w>>1] - 30));
            }

            else 
                for(int j = 0; j < kbottle; j++)
                    if(dst[w][j] < tmpd[j])
                    {
                        memcpy(dst[w], tmpd, sizeof(tmpd));
                        pq.push(nodeDistance(w, u, dst[w], cur.length + 1, cur.seq_len + unilen[w>>1] - 30));
                        break;
                    }
                    else if(dst[w][j] > tmpd[j])
                        break;

        }
    }

    leaves.clear();

    int sdis[kbottle];
    for(int i = 0; i < kbottle; i++)
        sdis[i] = INF;
    nodeDistance n(-1, -1, sdis, 0, 0);

    int x;
    vector<readinfo> ret;

    while(pq.size() > 0)
    {
        n = pq.top();
        pq.pop();
        
        x = n.id;

        while(vst[x] == start && !pq.empty())
        {
            n = pq.top();
            pq.pop();
            
            x = n.id;            
        }

        if(vst[x] == start && pq.empty())
            break;


        //if(start == 78)
            //printf("Dij %d %d %d %d %d\n",x, n.dis, n.pre, value[x], n.length);
        vst[x] = start;
        prevex[x] = n.pre;
        memcpy(dst[x], n.dis, sizeof(n.dis));


        if(n.seq_len > maxlen)
        {
            leaves.push_back(x);
            continue;
        }

        for (int i = h[x];i;i = q[i].next)
        {   
            int y = q[i].to;
            int tmpd[kbottle];

            for(int j = 0; j < kbottle; j++)
            {
                if(value[y] < dst[x][j])
                {
                    tmpd[j] = value[y];
                    for(int k = j + 1; k < kbottle; k++)
                        tmpd[k] = dst[x][k-1];
                    break;
                }

                tmpd[j] = dst[x][j];
            }

            if(vst[y] != start)
            {
                if(inpq[y] != start)
                {
                    inpq[y] = start;
                    memcpy(dst[y], tmpd, sizeof(tmpd));
                    pq.push(nodeDistance(y, x, dst[y], n.length + 1, n.seq_len + unilen[y>>1] - 30));
                }

                else
                {
                    for(int j = 0; j < kbottle; j++)
                        if(dst[y][j] < tmpd[j])
                        {
                            memcpy(dst[y], tmpd, sizeof(tmpd));
                            pq.push(nodeDistance(y, x, dst[y], n.length + 1, n.seq_len + unilen[y>>1] - 30));
                            break;
                        }
                        else if(dst[y][j] > tmpd[j])
                            break;
                }
            }
        }
    }

    for(int i = 0; i < v.size(); i++)
        if(vst[v[i].t] == start)
        {
            v[i].found = true;
            v[i].path = path2string(v[i].s, v[i].t, 31, v[i].l, v[i].r);
            if(v[i].path != "-1")
                ret.push_back(v[i]);
        }
    return ret;
}

int ReverseVertex(int x)
{
    if(x & 1)
        return x - 1;
    return x + 1;
}

vector<readinfo> all[maxn<<1];
int in[maxn<<1] = {0}, finish[maxn<<1] = {0};
queue<int> s_queue;

void Release(int x)
{
    finish[x] = 1;
    for(int i = h[x]; i; i = q[i].next)
        if(!finish[q[i].to])
        {
            in[q[i].to]--;
            if(in[q[i].to] == 0)
                s_queue.push(q[i].to);
        }
}

int Find_next(int s)
{
    queue<int> bfs_queue;

    for (int i = h[s]; i; i = q[i].next)
    {  
        int y = q[i].to;
        if(vst[y] == s && prevex[y] == s && !finish[y])
            bfs_queue.push(y);
    }



    while(!bfs_queue.empty())
    {
        int x = bfs_queue.front();
        bfs_queue.pop();

        if(all[x].size() > 0)
            return x;
        Release(x);

        for(int i = h[x]; i; i = q[i].next)
        {  
            int y = q[i].to;

            if(vst[y] == s && prevex[y] == x && !finish[y])
                bfs_queue.push(y);
        }
    }

    return -1;
}

int main(int argc, char *argv[])
{
    srand(time(NULL));

    ifstream input(argv[1]);
    string line;

    getline(input, line);
    stringstream fline(line);

    string tmp;
    getline(fline, tmp, '\t');
    n = stoi(tmp);
    getline(fline, tmp, '\t');
    m = stoi(tmp);
    int *hmap;
    hmap = new int[maxN];
    memset(hmap, 0, sizeof(hmap));

    for(int i = 0; i < n; i++)
    {
        getline(input, line);
        stringstream lines(line);

        string id, count, unitig;

        getline(lines, id, '\t');
        getline(lines, count, '\t');
        getline(lines, unitig, '\t');

        unis.push_back(unitig);
        value.push_back(stoi(count));
        unilen.push_back(unitig.length());
        int k = stoi(id);
        hmap[k] = i;
    }

    for(int i = 0; i < m; i++)
    {
        getline(input, line);
        stringstream lines(line);

        string a, b;

        getline(lines, a, '\t');
        getline(lines, b, '\t');
        
        int x = stoi(a);
        int y = stoi(b);

        x = (hmap[x>>1] << 1) + (x & 1);
        y = (hmap[y>>1] << 1) + (y & 1);
        Add(x,y);
    }
    
    int tmpcnt = 0;
    printf("Graph Build Finish Vertex = %d, Edge = %d\n", n<<1, edgecnt);
    
    ifstream query(argv[2]);

    for(int i = 0; i < (n<<1); i++)
    {
        vst[i] = -1;
        inpq[i] = -1;
    }

    int prevs = -1;
    vector<readinfo> v;
    vector<readinfo> ans;
    int qnum = 0, dnum = 0, pdnum = 0;
    int destnum = 0;

    while(1)
    {

        if(!getline(query, line))
            break;
        
        stringstream lines(line);

        string ids, ss, ts, l, r;
        int id, s, t;
        getline(lines, ids, '\t');
        getline(lines, ss, '\t');
        getline(lines, ts, '\t');
        getline(lines, l, '\t');
        getline(lines, r, '\t');
        qnum++;

        id = stoi(ids);
        s = stoi(ss);
        t = stoi(ts);

        s = (hmap[s>>1] << 1) + (s & 1);
        t = (hmap[t>>1] << 1) + (t & 1);

        readinfo tmp;
        tmp.s = s;
        tmp.t = t;
        tmp.id = id;
        tmp.l = l;
        tmp.r = r;
        prevs = s;
  
        all[s].push_back(tmp);
    }


    int undo = 0;

    for(int i = 0; i < (n<<1); i++)
        for(int j = h[i]; j ; j = q[j].next)
            in[q[j].to]++;

    for(int i = 0; i < (n<<1); i++)
    {
        if(in[i] == 0)
            s_queue.push(i);
        if(all[i].size() > 0)
            undo++;
    }

    qnum = 0;
    int pathnum = 0;
    clock_t start,end;
    start = clock();
    int c;
    int last = 0;
    printf("Query Nodes = %d\n", undo);
    while(undo)
    {
        int s = -1;
        if(!s_queue.empty())
        {
            s = s_queue.front();
            s_queue.pop();
            while(!s_queue.empty() && finish[s])
            {
                s = s_queue.front();
                s_queue.pop();
            }
        }

        if(s_queue.empty() && (s == -1 || finish[s]))
            for(int i = last; i < (n<<1); i++)
                if(!finish[i] && all[i].size())
                {
                    s = i;
                    last = i + 1;
                    break;
                }

        if(s == -1)
            break;

        if(all[s].size() == 0)
        {
            Release(s);
            continue;
        }

        v = Djikstra(s, all[s]);

        qnum += all[s].size();
        undo--;
        dnum++;

        for(int i = 0; i < v.size(); i++)
            if(v[i].found == 1)
                ans.push_back(v[i]);

        Release(s);

        while(h[s])
        {
            int nexts = Find_next(s);

            if(nexts == -1 || finish[nexts])
                break;

           
            v = PartlyDjikstra(nexts, s, all[nexts]);
            // Djikstra2(nexts, all[nexts], v.size());

            // for(int i = 0; i < (n<<1); i++)
            //     if((vst2[i] == nexts && nexts != vst[i])|| (vst[i] == nexts && nexts != vst2[i]) || (vst[i] == vst2[i] && vst[i] == nexts && dst[i] != dst2[i]))
            //     {
            //         printf("%d %d %d %d %d %d \n", nexts, i, dst[i], dst2[i], vst[i], vst2[i]);
            //     }

            // if(nexts == 78)
            //         cin>>c;

            qnum += all[nexts].size();
            if(all[nexts].size() > 0)
                undo--;
            pdnum++;

            for(int i = 0; i < v.size(); i++)
                if(v[i].found == 1)
                    ans.push_back(v[i]);
       

            s = nexts;
            Release(s);
        }

        pathnum++;
        if(pathnum % 10000 == 0)
            printf("Query: %d Ans: %d radio: %.3lf Nodes: %d Dij times: %d partlyDij times: %d Undo: %d\n", qnum, ans.size(), ans.size()/float(qnum), vnodes, dnum, pdnum, undo);
    }
    

    printf("Query: %d Ans: %d radio: %.3lf Nodes: %d Dij times: %d partlyDij times: %d \n", qnum, ans.size(), ans.size()/float(qnum), vnodes, dnum, pdnum);
    printf("Unfinshed %d\n", undo);
    end = clock();
    printf("Running time: %.2lf\n", (double)(end-start)/CLOCKS_PER_SEC/60);


    ofstream outq(argv[3]);
    for(int i = 0; i < ans.size(); i++)
        outq<<ans[i].id<<"\t"<<ans[i].path<<"\n";
    outq.close();

    return 0;
}