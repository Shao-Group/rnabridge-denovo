#include <iostream>
#include <stack>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "include/bifrost/CompactedDBG.hpp"
#include "include/bifrost/ColoredCDBG.hpp"

using namespace std;
const int INF = 0x7fffffff/2;
const int maxN = 20000000;
const int kmer = 31;
vector<string> dict;
int pokemon = 0;
string bridgeoutputpath;
string queryoutputpath;

void parse_ProgramOptions(int argc, char **argv, CCDBG_Build_opt& opt) {

    int option_index = 0, c;

    //const char* opt_string = "s:r:g:o:t:k:m:n:N:b:B:l:w:u:f:idvcya";
    const char* opt_string = "q:x:s:r:g:o:t:k:m:b:B:l:w:u:f:idvcya";

    static struct option long_options[] = {

        {"output-bridge",     required_argument,  0, 'q'},
        {"output-query",     required_argument,  0, 'x'},
        {"input-seq-files",     required_argument,  0, 's'},
        {"input-ref-files",     required_argument,  0, 'r'},
        {"input-graph-file",    required_argument,  0, 'g'},
        {"output-file",         required_argument,  0, 'o'},
        {"threads",             required_argument,  0, 't'},
        {"kmer-length",         required_argument,  0, 'k'},
        {"min-length",          required_argument,  0, 'm'},
        //{"num-kmers",           required_argument,  0, 'n'},
        //{"num-kmers2",          required_argument,  0, 'N'},
        {"bloom-bits",          required_argument,  0, 'b'},
        {"bloom-bits2",         required_argument,  0, 'B'},
        {"load-mbbf",           required_argument,  0, 'l'},
        {"write-mbbf",          required_argument,  0, 'w'},
        {"chunk-size",          required_argument,  0, 'u'},
        {"input-color-file",    required_argument,  0, 'f'},
        {"clip-tips",           no_argument,        0, 'i'},
        {"del-isolated",        no_argument,        0, 'd'},
        {"verbose",             no_argument,        0, 'v'},
        {"colors",              no_argument,        0, 'c'},
        {"keep-mercy",          no_argument,        0, 'y'},
        {"fasta",               no_argument,        0, 'a'},
        {0,                     0,                  0,  0 }
    };

    if (strcmp(argv[1], "build") == 0) opt.build = true;
    else if (strcmp(argv[1], "update") == 0) opt.update = true;

    if (opt.build || opt.update){

        while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {
            switch (c) {

                case 'q': {
                    bridgeoutputpath = optarg;
                    break;
                }

                case 'x': {
                    queryoutputpath = optarg;
                    break;
                }
                
                case 's': {

                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind){

                          opt.filename_seq_in.push_back(argv[optind]);
                    }

                    break;
                }
                case 'r': {

                    for (--optind; (optind < argc) && (*argv[optind] != '-'); ++optind){

                          opt.filename_ref_in.push_back(argv[optind]);
                    }

                    break;
                }
                case 'g':
                    opt.filename_graph_in = optarg;
                    break;
                case 'f':
                    opt.filename_colors_in = optarg;
                    break;
                case 'o':
                    opt.prefixFilenameOut = optarg;
                    break;
                case 't':
                    opt.nb_threads = atoi(optarg);
                    break;
                case 'k':
                    opt.k = atoi(optarg);
                    break;
                case 'm':
                    opt.g = atoi(optarg);
                    break;
                /*case 'n':
                    opt.nb_unique_kmers = atoi(optarg);
                    break;
                case 'N':
                    opt.nb_non_unique_kmers = atoi(optarg);
                    break;*/
                case 'b':
                    opt.nb_bits_unique_kmers_bf = atoi(optarg);
                    break;
                case 'B':
                    opt.nb_bits_non_unique_kmers_bf = atoi(optarg);
                    break;
                case 'w':
                    opt.outFilenameBBF = optarg;
                    break;
                case 'l':
                    opt.inFilenameBBF = optarg;
                    break;
                case 'u':
                    opt.read_chunksize = atoi(optarg);
                    break;
                case 'i':
                    opt.clipTips = true;
                    break;
                case 'd':
                    opt.deleteIsolated = true;
                    break;
                case 'v':
                    opt.verbose = true;
                    break;
                case 'c':
                    opt.outputColors = true;
                    break;
                case 'y':
                    opt.useMercyKmers = true;
                    break;
                case 'a':
                    opt.outputGFA = false;
                    break;
                default: break;
            }
        }
    }
}

bool check_ProgramOptions(CCDBG_Build_opt& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    auto check_files = [&](vector<string>& v_files) {

        vector<string> files_tmp;

        char* buffer = new char[4096]();

        for (const auto& file : v_files) {

            if (!check_file_exists(file)) {

                cerr << "Error: File " << file << " not found." << endl;
                ret = false;
            }
            else {

                const string s_ext = file.substr(file.find_last_of(".") + 1);

                if ((s_ext == "txt")){

                    FILE* fp = fopen(file.c_str(), "r");

                    if (fp != NULL){

                        fclose(fp);

                        ifstream ifs_file_txt(file);
                        istream i_file_txt(ifs_file_txt.rdbuf());

                        while (i_file_txt.getline(buffer, 4096)){

                            fp = fopen(buffer, "r");

                            if (fp == NULL) {

                                cerr << "Error: Could not open file " << buffer << " for reading." << endl;
                                ret = false;
                            }
                            else {

                                fclose(fp);
                                files_tmp.push_back(string(buffer));
                            }
                        }

                        ifs_file_txt.close();
                    }
                    else {

                        cerr << "Error: Could not open file " << file << " for reading." << endl;
                        ret = false;
                    }
                }
                else files_tmp.push_back(file);
            }
        }

        v_files = move(files_tmp);

        delete[] buffer;
    };

    // Check general parameters

    if (!opt.build && !opt.update){

        cerr << "Error: No command selected (can be 'build' or 'update')." << endl;
        ret = false;
    }

    if (opt.nb_threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.nb_threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << "." << endl;
        ret = false;
    }

    if (opt.k <= 0){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << "." << endl;
        ret = false;
    }

    if (opt.g <= 0){

        cerr << "Error: Length m of minimizers cannot be less than or equal to 0." << endl;
        ret = false;
    }

    if (opt.g > opt.k - 2){

        cerr << "Error: Length m of minimizers cannot exceed k - 2 (" << (opt.k - 2) << ")." << endl;
        ret = false;
    }

    const string out = opt.prefixFilenameOut + (opt.outputGFA ? ".gfa" : ".fasta");

    FILE* fp = fopen(out.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << out << "." << endl;
        ret = false;
    }
    else {

        fclose(fp);
        if (remove(out.c_str()) != 0) cerr << "Error: Could not remove temporary file " << out << "." << endl;
    }

    if ((opt.filename_seq_in.size() + opt.filename_ref_in.size()) == 0) {

        cerr << "Error: Missing input files." << endl;
        ret = false;
    }
    else {

        check_files(opt.filename_seq_in);
        check_files(opt.filename_ref_in);
    }

    if (opt.build){ // Check param. command build

        if (opt.read_chunksize <= 0) {

            cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0." << endl;
            ret = false;
        }

        if (opt.outFilenameBBF.length() != 0){

            FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

            if (fp == NULL) {

                cerr << "Error: Could not open Blocked Bloom filter file " << opt.outFilenameBBF << " for writing." << endl;
                ret = false;
            }
            else {

                fclose(fp);

                if (remove(opt.outFilenameBBF.c_str()) != 0){

                    cerr << "Error: Could not remove temporary file " << opt.outFilenameBBF << "." << endl;
                }
            }
        }

        if (opt.inFilenameBBF.length() != 0){

            if (check_file_exists(opt.inFilenameBBF)){

                FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input Blocked Bloom filter file " << opt.inFilenameBBF << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
            else {

                cerr << "Error: Input Blocked Bloom filter " << opt.inFilenameBBF << " file does not exist." << endl;
                ret = false;
            }
        }
    }

    if (opt.update){

        if (opt.filename_graph_in.length() == 0){

            cerr << "Error: No graph file to update was provided in input." << endl;
            ret = false;
        }
        else if (!check_file_exists(opt.filename_graph_in)){

            cerr << "Error: The graph file to update does not exist." << endl;
            ret = false;
        }
        else {

            FILE* fp = fopen(opt.filename_graph_in.c_str(), "r");

            if (fp == NULL) {

                cerr << "Error: Could not read input graph file " << opt.filename_graph_in << "." << endl;
                ret = false;
            }
            else fclose(fp);
        }

        if (opt.filename_colors_in.length() != 0){

            if (!check_file_exists(opt.filename_colors_in)){

                cerr << "Error: The input color file does not exist." << endl;
                ret = false;
            }
            else {

                FILE* fp = fopen(opt.filename_colors_in.c_str(), "rb");

                if (fp == NULL) {

                    cerr << "Error: Could not read input color file " << opt.filename_colors_in << "." << endl;
                    ret = false;
                }
                else fclose(fp);
            }
        }
    }

    return ret;
}

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

        int subgraph = 0;
        int id = 0;
        int sub_id = 0;
    private:
        bool v;
};

void clearMarking(const unordered_set<UnitigColorMap<myUnitigData>, UnitigMapHash<DataAccessor<myUnitigData>, DataStorage<myUnitigData>>>& set_km_seen){

    for (const auto& ucm : set_km_seen){

        DataAccessor<myUnitigData>* da_ucm = ucm.getData();
        myUnitigData* data_ucm = da_ucm->getData(ucm);

        data_ucm->clear(ucm);
    }
}

void clearMarking(ColoredCDBG<myUnitigData>& ccdbg){

    for (const auto& unitig : ccdbg){

        DataAccessor<myUnitigData>* da_ucm = unitig.getData();
        myUnitigData* data_ucm = da_ucm->getData(unitig);

        data_ucm->clear(unitig);
    }
}

int kk = 0, global_id = 0;
int g_id = 0;

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
        if(a[i] == 'N')
            b = b + 'N';
    }

    return b;
}

int h[maxN<<1] = {0}, cnt = 0;
struct edgelist{int to,next;}
q[maxN<<2];

void Add(int a, int b)
{
    q[++cnt] = (edgelist){b, h[a]};
    h[a] = cnt;
}

bool Check_Direction(string a, string b)
{
    string a0 = a.substr(a.length() - 30, 30);
    string b0 = b.substr(0, 30);

    return a0 == b0;
}

int tp1n = 0, tp2n = 0, tp3n = 0, tp4n = 0;
int pn = 0;

void DFS_Iterative(const UnitigColorMap<myUnitigData>& ucm)
{
    stack<UnitigColorMap<myUnitigData>> stck;
    UnitigColorMap<myUnitigData> ucm_tmp(ucm);
    stck.push(ucm_tmp);
    
    DataAccessor<myUnitigData>* da = ucm_tmp.getData();
    myUnitigData* data = da->getData(ucm_tmp);
    
    dict.push_back("");
    data->set_visited();
    data->id = ++global_id;

    while (!stck.empty())
    { 
        ucm_tmp = stck.top();
        stck.pop();
        pn++;

        da = ucm_tmp.getData();
        data = da->getData(ucm_tmp);

        data->sub_id = ++kk;
        data->subgraph = g_id;
        dict[data->id] = ucm_tmp.referenceUnitigToString();
        //cout<<data->id<<" "<<dict[data->id]<<endl;

        string s1 = ucm_tmp.referenceUnitigToString();

        for (auto& successor : ucm_tmp.getSuccessors()) 
        {
            DataAccessor<myUnitigData>* su_da = successor.getData();
            myUnitigData* su_data = su_da->getData(successor);
           
            if(su_data->is_not_visited())
            {
                stck.push(successor);
                su_data->set_visited();
                su_data->id = ++global_id;
                dict.push_back("");
            }
            //else 
            //    continue;


            string s2 = successor.referenceUnitigToString();
            int id1 = data->id;
            int id2 = su_data->id;
            //cout<<"su "<<id2<<endl;

            if(Check_Direction(s1, s2))
            {
                //tp1n++;
                Add((id1 << 1) , (id2 << 1));
            }
            
            else if(Check_Direction(s1, trans(s2)))
            {
                //tp2n++;
                Add((id1 << 1), (id2 << 1) + 1);
            }

            else if(Check_Direction(s2, s1))
            {
                //tp3n++;
                Add((id1 << 1) + 1, (id2 << 1) + 1);
            }
            
            else if(Check_Direction(trans(s2), s1))
            {
                //tp4n++;
                Add((id1 << 1) + 1, (id2 << 1));
            }

            // else
            //     printf("%s\n%s\n%s\n%s\n\n",s1.c_str(), trans(s1).c_str(), s2.c_str(), trans(s2).c_str());
        }

        for (auto& predecessor : ucm_tmp.getPredecessors()) 
        {
            DataAccessor<myUnitigData>* su_da = predecessor.getData();
            myUnitigData* su_data = su_da->getData(predecessor);
           
            if(su_data->is_not_visited())
            {
                stck.push(predecessor);
                su_data->set_visited();
                su_data->id = ++global_id;
                dict.push_back("");
            }
            //else 
              //  continue;

            string s2 = predecessor.referenceUnitigToString();

            int id1 = data->id;
            int id2 = su_data->id;

            if(Check_Direction(s2, s1))
            {
                tp1n++;
                Add((id1 << 1) + 1, (id2 << 1) + 1);
            }
            
            else if(Check_Direction(s2, trans(s1)))
            {
                tp2n++;
                Add((id1 << 1), (id2 << 1) + 1);
            }

            else if(Check_Direction(s1, s2))
            {
                tp3n++;
                Add((id1 << 1), (id2 << 1));
            }
            
            else if(Check_Direction(trans(s1), s2))
            {
                tp4n++;
                Add((id1 << 1) + 1, (id2 << 1));
            }
        }
    }
}

size_t getNbConnectedComponent(ColoredCDBG<myUnitigData>& ccdbg)
{

    size_t nb_cc = 0; // Number of connected components
    dict.push_back("");

    for (auto& unitig : ccdbg){ // Iterate over unitigs of a colored de Bruijn graph

        DataAccessor<myUnitigData>* da = unitig.getData(); // Get DataAccessor from unitig
        myUnitigData* data = da->getData(unitig); // Get boolean from DataAccessor


        if (data->is_not_visited()){ // If boolean indicates the unitig was not visited
            ++nb_cc; // It's a new connected components
            ++g_id;
            kk = 0;

            DFS_Iterative(unitig);
        }
    }
    clearMarking(ccdbg);
    return nb_cc;
}


ColoredCDBG<myUnitigData> cdbg1;
vector<int> countv[maxN];
int val[maxN<<1] = {0};

void Init(string read, ColoredCDBG<myUnitigData>& ccdbg)
{
    int pos, length;
    

    length = read.length();
    pos = 0;
    const char* s = read.c_str();

    while(pos + kmer - 1 < length)
    {
        UnitigMap<DataAccessor<myUnitigData>, DataStorage<myUnitigData>> tmp = ccdbg.findUnitig(s, pos, length);
        if(tmp.isEmpty)
        {
            pos++;
        }
        else
        {
            int l = tmp.len;

            DataAccessor<myUnitigData>* sda = tmp.getData();
            myUnitigData* sdata = sda->getData(tmp);

            int k = sdata->id;
            if(countv[k].empty())
                for(int i = 0; i < tmp.size - kmer + 1; i++)
                    countv[k].push_back(0);

            for(int i = 0; i < l; i++)
                countv[k][tmp.dist + i]++;

            // string ref = tmp.referenceUnitigToString();
            // if(!tmp.strand)
            //     ref = trans(ref);
            // printf("%s\n%s\n%d %d\n", read.substr(pos, l + kmer - 1).c_str(), ref.c_str(), tmp.dist, l);
            pos = pos + l;

        }
    }
}



struct Eclass{
    string l,r;
    int lid, rid, id;
    //string unitigl, unitigr;
    Eclass* next;
};

Eclass* E_pointer[maxN<<1 + 5];   

void get_Eclass(string l, string r, int id, ColoredCDBG<myUnitigData>& ccdbg)
{
    int posh, post, length;
    int k_length = ccdbg.getK();

    length = l.length();
    posh = 0;
    post = length - k_length;

    const char* s = l.c_str();
    string refl, refr;
    //cout<<l<<" "<<r<<endl;
    UnitigMap<DataAccessor<myUnitigData>, DataStorage<myUnitigData>> lnode = ccdbg.findUnitig(s, post, length);

    if(lnode.isEmpty)
        return;
    const char* s1 = r.c_str();

    UnitigMap<DataAccessor<myUnitigData>, DataStorage<myUnitigData>> rnode = ccdbg.findUnitig(s1, posh, k_length);

    if(rnode.isEmpty)
        return;

    DataAccessor<myUnitigData>* lda = lnode.getData();
    myUnitigData* ldata = lda->getData(lnode);

    DataAccessor<myUnitigData>* rda = rnode.getData();
    myUnitigData* rdata = rda->getData(rnode);


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

string upper(string s)
{
    int len = s.length();
    string ret;
    for(int i = 0; i < len; i++)
        if(s[i] >= 'a' && s[i] <= 't')
            ret += s[i] - 'a' + 'A';
        else if(s[i] >= 'A' && s[i] <= 'T')
            ret += s[i];

    return ret;
}

int threshold = 20;
int main(int argc, char **argv)
{

    CCDBG_Build_opt opt;

    opt.outputColors = false; // We dont know yet if we want colors or not

    parse_ProgramOptions(argc, argv, opt); // Parse input parameters

    ColoredCDBG<myUnitigData> cdbg1(opt.k, opt.g);
    cdbg1.read(opt.filename_graph_in, opt.filename_colors_in, opt.nb_threads, opt.verbose);

    const size_t cdbg1_len = cdbg1.length();
    printf("Origin map size: %d\n",cdbg1.size());


    cout << "=== Computing number of connected components ===" << endl;

    const size_t nb_cc = getNbConnectedComponent(cdbg1);
    cout << nb_cc << " connected components found" << endl;
    vector<string> mycolors = cdbg1.getColorNames();

    printf("%d\n", mycolors.size());
    for(int i = 0; i < mycolors.size(); i++)
        printf("%s\n", mycolors[i].c_str());

    int k_length = cdbg1.getK();

    ifstream infile(mycolors[0].c_str());
    ifstream infile1(mycolors[1].c_str());

    string line, line1;
    int readsnb = 0;

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

        line = upper(line);
        line1 = upper(line1);
        
        Init(line, cdbg1);
        Init(line1, cdbg1);
        get_Eclass(line, trans(line1), readsnb, cdbg1);
        //get_Eclass(line1, trans(line), readsnb, cdbg1);

        getline(infile, line);
        getline(infile, line);
        getline(infile1, line1);
        getline(infile1, line1);
        if(readsnb % 10000000 == 0)
            printf("%d\n", readsnb);
    }
    printf("Initilize complete \n");


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

    for(int i = 2; i <= global_id * 2 + 1; i++)
        if(val[i] > threshold)
            for(int j = h[i]; j; j = q[j].next)
                if(val[q[j].to] >= threshold)
                    enumber++;

    ofstream outbri(bridgeoutputpath);
    outbri<<vnumber<<"\t"<<enumber<<endl;

    for(int i = 1; i <= global_id; i++)
        if(val[i<<1] >= threshold)
            outbri<<i<<"\t"<<val[i<<1]<<"\t"<<dict[i]<<endl;

    for(int i = 2; i <= global_id * 2 + 1; i++)
        if(val[i] > threshold)
            for(int j = h[i]; j; j = q[j].next)
                if(val[q[j].to] >= threshold)
                    outbri<<i<<"\t"<<q[j].to<<"\n";
                

    outbri.close();
    
    ofstream outq(queryoutputpath);

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
    return 0;      
}