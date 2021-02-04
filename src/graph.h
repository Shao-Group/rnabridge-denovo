#include "bifrost/CompactedDBG.hpp"
#include "bifrost/ColoredCDBG.hpp"
#include "utils.h"

using namespace std;

#ifndef __GRAPH_H__
#define __GRAPH_H__

class myUnitigData : public CCDBG_Data_t<myUnitigData>, CDBG_Data_t<myUnitigData> {

    public:

        myUnitigData() : v(0){} // Initiate the boolean to "not visited"

        // Clear method for CompactedDBG
        void clear(const UnitigMap<myUnitigData>& um_dest)
        {
            set_not_visited(); // Set the new unitig to "not visited"
        }

         // Clear method for CompactedDBG
        void clear(const UnitigColorMap<myUnitigData>& um_dest)
        {
            set_not_visited(); // Set the new unitig to "not visited"
        }

        inline void set_visited() { v = 1; }
        inline void set_not_visited() { v = 0; }

        inline bool is_visited() const { return (v == 1); }
        inline bool is_not_visited() const { return (v == 0); }

        int id = 0;
    private:
        bool v;
};


class graph{

    const int INF = 0x7fffffff/2;
    static const int maxN = 20000000;
    int kmer = 31;
    int threshold = 1;
    int global_id = 0;

    vector<string> dict;
    vector<int> countv[maxN];
    int val[maxN<<1] = {0};

    int h[maxN<<1] = {0}, cnt = 0;
    struct edgelist{int to,next;} q[maxN<<2];

    struct Eclass{
        string l,r;
        int lid, rid, id;
        //string unitigl, unitigr;
        Eclass* next;
    };
    Eclass* E_pointer[maxN<<1 + 5];   


    
    CCDBG_Build_opt opt, opt1;

    ColoredCDBG<> dbg;
    ColoredCDBG<myUnitigData> cdbg;

    void init(string read, ColoredCDBG<myUnitigData>& ccdbg);
    void get_Eclass(string l, string r, int id, ColoredCDBG<myUnitigData>& ccdbg);
    void add(int a, int b);
    void dfs_Iterative(const UnitigColorMap<myUnitigData>& ucm);
    void dfs();
public:
    void build(string rd1, string rd2, string output);

};

#endif