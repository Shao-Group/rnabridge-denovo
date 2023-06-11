#include "graph.h"
#include <stack>

bool check_Direction(string a, string b)
{
    int len = a.length();

    for(int i = 0; i < 30; i++)
        if(b[i] != a[len-30+i])
            return 0;

    return 1;
}


void upper(string &s)
{
    int len = s.length();
    for(int i = 0; i < len; i++)
        if(s[i] >= 'a' && s[i] <= 't')
            s[i] += 'A' - 'a';
}

void graph::build(string rd1, string rd2, string output)
{
    printf("Read de Bruijn Graph\n");
    // opt.clipTips = false;
    // opt.deleteIsolated = false;
    // opt.filename_seq_in.push_back(rd1);
    // opt.filename_seq_in.push_back(rd2);
    // opt.nb_threads = 15;
    // opt.build = true;
    // opt.prefixFilenameOut = output + "/dbg";

    // if(!dbg.build(opt))
    //     printf("Build de Bruijn Graph Failed with Bifrost\n");
    
    // if(!dbg.write(opt.prefixFilenameOut, 15))
    //     printf("Fail to write graph to %s\n", opt.prefixFilenameOut.c_str());


    opt1.nb_threads = 15;
    opt1.update = true;
    opt1.filename_graph_in = output + "/dbg.gfa";
    cdbg.read(opt1.filename_graph_in, opt1.nb_threads, opt1.verbose);

    printf("Origin graph size: %d\n", int(cdbg.size()));
    dfs();
    cout<<edges<<endl;
    
    string line, line1;
    int readsnb = 0;
    ifstream infile(rd1.c_str());
    ifstream infile1(rd2.c_str());

    while(1)
    {
        if(!getline(infile, line))
            break;

        if(!getline(infile, line))
            break;

        if(!getline(infile1, line1))
            break;

        if(!getline(infile1, line1))
            break;

        readsnb++;
        //cout<<readsnb<<" "<<line<<" "<<line1<<endl;
        upper(line);
        upper(line1);
        
        init(line, cdbg);
        init(line1, cdbg);
        get_Eclass(line, trans(line1), readsnb, cdbg);

        getline(infile, line);
        getline(infile, line);
        getline(infile1, line1);
        getline(infile1, line1);
        if(readsnb % 10000000 == 0)
            printf("Aligning reads: %d\n", readsnb);
    }
    printf("Aligning reads to graph complete \n");


    for(int i = 1; i <= global_id; i++)
    {
        int minv = 0x7fffffff/2;
        for(auto j: countv[i])
            if(minv > j)
                minv = j;
        if(minv == 0x7fffffff/2)
            minv = 0;

        val[i<<1] = val[(i<<1) + 1] = minv;
    }
    int vnumber = 0, enumber = 0;

    for(int i = 1; i <= global_id; i++)
        if(val[i<<1] >= threshold)
            vnumber++;

    printf("Vertices: %d\n", vnumber);
    for(int i = 2; i <= global_id * 2 + 1; i++)
        if(val[i] >= threshold)
            for(int j = h[i]; j; j = q[j].next)
                if(val[q[j].to] >= threshold)
                    enumber++;

    printf("Edges: %d\n", enumber);

    ofstream outbri(output + "/graph");

    outbri<<vnumber<<"\t"<<enumber<<endl;

    for(int i = 1; i <= global_id; i++)
        if(val[i<<1] >= threshold)
            outbri<<i<<"\t"<<val[i<<1]<<"\t"<<dict[i]<<endl;

    for(int i = 2; i <= global_id * 2 + 1; i++)
        if(val[i] >= threshold)
            for(int j = h[i]; j; j = q[j].next)
                if(val[q[j].to] >= threshold)
                    outbri<<i<<"\t"<<q[j].to<<"\n";
                

    outbri.close();

    ofstream outq(output + "/query");

    for(int i = 2; i <= global_id * 2 + 1; i++)
        if(E_pointer[i] != nullptr && val[i] >= threshold)
        {       

            Eclass* e = E_pointer[i];

            while(e != nullptr)
            {
                if(val[e->rid] >= threshold)
                    outq<<e->id<<"\t"<<e->lid<<"\t"<<e->rid<<"\t"<<e->l<<"\t"<<e->r<<endl;
                e = e->next;
            }
    }

    outq.close();
    printf("Weighted graph saved\n");

}

void graph::dfs()
{
    dict.push_back("");

    for (auto& unitig : cdbg)// Iterate over unitigs of a colored de Bruijn graph
    { 

        myUnitigData* data = unitig.getData(); // Get DataAccessor from unitig

        if(data->is_not_visited())// If boolean indicates the unitig was not visited
            dfs_Iterative(unitig);
    }
}



void graph::dfs_Iterative(const UnitigMap<myUnitigData>& ucm)
{
    stack<UnitigMap<myUnitigData>> stck;
    UnitigMap<myUnitigData> ucm_tmp(ucm);
    stck.push(ucm_tmp);
    
   
    myUnitigData* data = ucm_tmp.getData();
    
    dict.push_back(ucm_tmp.referenceUnitigToString());
    data->set_visited();
    data->id = ++global_id;

    while (!stck.empty())
    { 
        ucm_tmp = stck.top();
        stck.pop();

        data = ucm_tmp.getData();
        string s1 = dict[data->id];

        for (auto& successor : ucm_tmp.getSuccessors()) 
        {
            myUnitigData* su_data = successor.getData();
            edges++;

            if(su_data->is_not_visited())
            {
                stck.push(successor);
                su_data->set_visited();
                su_data->id = ++global_id;
                dict.push_back(successor.referenceUnitigToString());
            }


            string s2 = successor.referenceUnitigToString();

            int id1 = data->id;
            int id2 = su_data->id;

            if(check_Direction(s1, s2))
                add((id1 << 1) , (id2 << 1));
            
            else if(check_Direction(s1, trans(s2)))
                add((id1 << 1), (id2 << 1) + 1);

            else if(check_Direction(s2, s1))
                add((id1 << 1) + 1, (id2 << 1) + 1);

            else if(check_Direction(trans(s2), s1))
                add((id1 << 1) + 1, (id2 << 1));
        }

        for (auto& predecessor : ucm_tmp.getPredecessors()) 
        {
            myUnitigData* su_data = predecessor.getData();
            edges++;

            if(su_data->is_not_visited())
            {
                stck.push(predecessor);
                su_data->set_visited();
                su_data->id = ++global_id;
                dict.push_back(predecessor.referenceUnitigToString());
            }

            string s2 = predecessor.referenceUnitigToString();

            int id1 = data->id;
            int id2 = su_data->id;

            if(check_Direction(s2, s1))
                add((id1 << 1) + 1, (id2 << 1) + 1);

            else if(check_Direction(s2, trans(s1)))
                add((id1 << 1), (id2 << 1) + 1);

            else if(check_Direction(s1, s2))
                add((id1 << 1), (id2 << 1));
            
            else if(check_Direction(trans(s1), s2))
                add((id1 << 1) + 1, (id2 << 1));
        }
    }
}


void graph::add(int a, int b)
{
    q[++cnt] = (edgelist){b, h[a]};
    h[a] = cnt;
}


void graph::init(string read, CompactedDBG<myUnitigData>& ccdbg)
{
    int pos, length;
    

    length = read.length();
    pos = 0;
    const char* s = read.c_str();

    while(pos + kmer - 1 < length)
    {
        auto tmp = ccdbg.findUnitig(s, pos, length);
        if(tmp.isEmpty)
        {
            pos++;
        }
        else
        {
            int l = tmp.len;

            myUnitigData* sdata = tmp.getData();

            int k = sdata->id;
            if(countv[k].empty())
                for(int i = 0; i < int(tmp.size) - kmer + 1; i++)
                    countv[k].push_back(0);

            for(int i = 0; i < l; i++)
                countv[k][tmp.dist + i]++;

            pos = pos + l;
        }
    }
}

void graph::get_Eclass(string l, string r, int id, CompactedDBG<myUnitigData>& ccdbg)
{
    int posh, post, length;
    int k_length = ccdbg.getK();

    length = l.length();
    posh = 0;
    post = length - k_length;

    const char* s = l.c_str();
    string refl, refr;
    //cout<<l<<" "<<r<<endl;
    auto lnode = ccdbg.findUnitig(s, post, length);

    if(lnode.isEmpty)
        return;
    const char* s1 = r.c_str();

    auto rnode = ccdbg.findUnitig(s1, posh, k_length);

    if(rnode.isEmpty)
        return;

    myUnitigData* ldata = lnode.getData();

    myUnitigData* rdata = rnode.getData();


    string stringl = lnode.referenceUnitigToString();
    string stringr = rnode.referenceUnitigToString();
    int start, destination;

    if(stringl.find(l.substr(post, k_length)) == string::npos)
        start = (ldata->id << 1) + 1;
    else
        start = ldata->id << 1;

    if(stringr.find(r.substr(posh, k_length)) == string::npos)
        destination = (rdata->id << 1) + 1;
    else
        destination  = rdata->id << 1;


    Eclass* k = new Eclass;

    if(E_pointer[start] != nullptr)
        k->next = E_pointer[start];
    else
        k->next = nullptr;

    E_pointer[start] = k;

    k->lid = start;
    k->rid = destination;
    k->id = id;
    k->l = l;
    k->r = r;
}

