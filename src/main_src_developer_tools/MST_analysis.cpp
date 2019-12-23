/**********************************************************
  *
  *
  *
  *
  *
  ********************************************************/
#include "algorithm/disjoin_set/disjoin_set.h"

#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "stLFR/CBB.h"

#include <iostream>
#include <sstream>
#include <regex>

#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <tuple>

/********************************************************
  *
  * basic structures
  *
  ******************************************************/

typedef BGIQD::Algorithm::DisJoin_Set<unsigned int> DisJoin_Set;
DisJoin_Set dis_joint_set;

struct MST_AnalysisEdge {
    unsigned int from ;
    unsigned int to ;
    int rank ; // 0 means non-unique relative edge
    float sim ;
};

struct MST_AnalysisNode ;

MST_AnalysisNode & get_oppo( int edge_id , unsigned int node_id ) ;

struct MST_AnalysisNode {

    unsigned int id ;
    bool operator < ( const MST_AnalysisNode & i ) const { return id < i.id ; } 
    bool is_unique ;

    enum GraphType {
        UnknowNode = 0 ,
        LinearNode = 1 ,
        TipNode = 2 ,
        JunctionNode = 3 ,
    };
    GraphType graph_type ;
    enum JunctionType {
        UnknowJunction = 0 ,
        TipJunction = 1,
        LongJunction = 2
        //BrankJunction = 2 ,
        //MixedJunction = 3 ,
        //NeibJunction = 4 
    };

    enum NodeL2 {
        UnknowL2 = 0 ,
        LowL2 = 1 ,
        HighL2 = 2 
    };

    enum NeibType { 
        UnknowNieb = 0 ,
        SimpleNieb = 1 ,
        Short_Long = 2 ,
        Long_Log   = 3 ,
        MaxNieb    = 4
    };

    JunctionType junction_type ;

    NeibType nieb_type ;

    NodeL2 l2_type ;

    std::vector<int> edges ; // only store the index of edges .

    int BranchNum () const { return edges.size() ; }

    void AddEdge(unsigned int self_id , int index)
    {
        id = self_id ;
        is_unique = false ;
        graph_type = UnknowNode ;
        junction_type = UnknowJunction ;
        edges.push_back(index);
    }

    std::string ref ;
    int pos;
    int rank ;

            // id <-->      len
    std::map<unsigned int , int >  branch_len_map ;
    int max_branch_len ;
    unsigned int max_branch_from_id ;

    std::set<unsigned int> SendDepth() {
        std::set<unsigned int> ret ;
        for( auto x : edges ) {
            auto & nieb_node = get_oppo(x,id);
            if( nieb_node.id != max_branch_from_id ) 
                if( nieb_node.UpdateMaxBranchNode(max_branch_len,id) )
                    ret.insert(nieb_node.id);
        }
        return ret ;
    }
    bool UpdateMaxBranchNode(int len , unsigned int id ) {
        branch_len_map[id] = len+1 ;
        if( len+1 > max_branch_len ) {
            max_branch_len = len +1;
            max_branch_from_id = id ;
            return true ;
        }
        return false ;
    }
};

struct QuastInfo
{
    int ref_start ;
    int ref_end ;
    int query_start ;
    int query_end ;
    std::string ref ;
    unsigned int contig ;

    void InitFromString(const std::string & line )
    {
        std::string tmp1 ;
        std::istringstream ist(line);
        ist>>ref_start>>ref_end>>query_start>>query_end>>ref>>tmp1;
        std::string ctg ;
        for( auto i : tmp1 )
        {
            if( std::isdigit(i) )
            {
                ctg += i;
            }
            else
                break ;
        }
        contig = std::stoul(ctg);
    }
};

/********************************************************
  *
  * global variables
  *
  ******************************************************/

static std::string output_dir;
std::map<unsigned int , MST_AnalysisNode> nodes;
std::vector<MST_AnalysisEdge> edges;
std::map<unsigned int , std::map<unsigned int , float > > JS_matrix;
//       ref                       pos   contig
std::map<std::string , std::vector<std::tuple<int , unsigned int> > > ref_ranks ;
/********************************************************
  *
  * logic functions
  *
  ******************************************************/

int parse_mintree( std::istream & ist )
{
    int ret = 0;
    std::string line ;
    while( ! ist.eof() )
    {
        std::getline(ist,line);
        if(line.empty() ) continue ;
        if(line == "graph {") continue ;
        if(line == "}") continue ;
        unsigned long long from , to ;
        float sim ;
        std::istringstream isst(line);
        //17      --      2849 [ label="0.367841" ]
        sscanf(line.c_str(), "%llu\t--\t%llu [ label=\"%f\" ]",&from, &to,&sim);
        edges.push_back({ (unsigned int )from , (unsigned int)to , 0 , sim});
        nodes[from].AddEdge(from,edges.size()-1);
        nodes[to].AddEdge(to ,edges.size()-1);
        dis_joint_set.AddConnect(from , to );
        ret ++ ;
    }
    return ret ;
}

int parse_sort_unique_contigs(std::istream & ist)
{
    std::string line ;
    int ret = 0 ;
    while( ! ist.eof() )
    {
        std::getline(ist,line);
        if(line.empty() ) continue ;
        QuastInfo tmp ;
        tmp.InitFromString(line);
        if( nodes.find(tmp.contig) != nodes.end() )
        {
            nodes[tmp.contig].is_unique = true ;
            nodes[tmp.contig].ref = tmp.ref ;
            nodes[tmp.contig].pos = tmp.ref_start ;
            ret ++ ;
        }
    }
    return ret;
}

MST_AnalysisNode & get_oppo( int edge_index , unsigned int self_id )
{
    auto & the_edge = edges.at(edge_index) ;
    if( the_edge.from == self_id ) 
        return nodes.at(the_edge.to);
    else if ( the_edge.to == self_id )
        return nodes.at(the_edge.from);
    else
    {
        assert(0);
        static MST_AnalysisNode err;
        return err;
    }
}
int  get_rank( unsigned int from , unsigned int to) 
{
    if(nodes.find(from) == nodes.end()) return 0 ;
    if(nodes.find(to) == nodes.end() ) return 0 ;
    const auto & from_node = nodes.at(from);
    const auto & to_node = nodes.at(to);
    if( ! from_node.is_unique  || from_node.ref == "" ) return 0 ;
    if( ! to_node.is_unique  || to_node.ref == "" ) return 0 ;
    if( from_node.ref != to_node.ref )  return -2 ;
    return std::abs(from_node.rank  - to_node.rank );
}

int count_edge_rank(int num)
{
    std::string file ;
    if( num > 0 )
        file = output_dir+"/edge_rank_"+std::to_string(num)+".txt";
    else if (num == 0 )
        file = output_dir+"/edge_non_unique.txt";
    else if (num == -2 )
        file = output_dir+"/edge_cross_ref.txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    int ret = 0 ;
    for( const auto & edge : edges )
        if(edge.rank == num ) 
        {
            ret ++ ;
            (*out)<<edge.from<<'\t'<<edge.to<<'\t'<<edge.sim<<'\n';
        }
    delete out ;
    return ret ;
}
int count_edge_rank_gt(int num)
{
    int ret = 0 ;
    std::string file = output_dir+"/edge_rank_gt_"+std::to_string(num)+".txt";
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto & edge : edges )
        if(edge.rank > num ) 
        {
            ret ++ ;
            (*out)<<edge.from<<'\t'<<edge.to<<'\t'<<edge.sim<<'\n';
        }
    delete out ;
    return ret ;
}

int count_linear()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::LinearNode) ret ++ ;
    }
    return ret ;
}
int count_linear_edge(const MST_AnalysisNode & node )
{
    int ret = 0 ;
    for( int index : node.edges )
        if( edges.at(index).rank == 1 )
            ret ++ ;
    return ret ;
}
int count_linear_edge_rank1(int num)
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::LinearNode
           && count_linear_edge(i) == num )
            ret ++ ;
    }
    return ret ;
}

int count_tip()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::TipNode) ret ++ ;
    }
    return ret ;
}
int count_junction_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::JunctionNode  && ! i.is_unique ) ret ++ ;
    }
    return ret ;
}
int count_junction_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::JunctionNode  && i.is_unique ) ret ++ ;
    }
    return ret ;
}
int count_junction()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.graph_type == MST_AnalysisNode::GraphType::JunctionNode ) ret ++ ;
    }
    return ret ;
}
int count_junction(int junc_num)
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.BranchNum() == junc_num  ) ret ++ ;
    }
    return ret ;
}

int count_junction_gt(int junc_num)
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.BranchNum() > junc_num  ) ret ++ ;
    }
    return ret ;
}

int count_long_junction()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::LongJunction ) ret ++ ;
    }
    return ret ;
}

int count_long_junction_low()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::LongJunction
            && i.l2_type == MST_AnalysisNode::NodeL2::LowL2 
                ) ret ++ ;
    }
    return ret ;
}

int count_long_junction_high()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::LongJunction
            && i.l2_type == MST_AnalysisNode::NodeL2::HighL2
                ) ret ++ ;
    }
    return ret ;
}

int count_tip_junction() 
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::TipJunction ) ret ++ ;
    }
    return ret ;
}

int count_tip_junction_low()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::TipJunction
            && i.l2_type == MST_AnalysisNode::NodeL2::LowL2 
                ) ret ++ ;
    }
    return ret ;
}
int count_tip_junction_high()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::TipJunction
            && i.l2_type == MST_AnalysisNode::NodeL2::HighL2
                ) ret ++ ;
    }
    return ret ;
}
int count_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( ! i.is_unique ) ret ++ ;
    }
    return ret ;
}

void prepare_ref_ranks()
{
    // collect ref data
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.is_unique ) ref_ranks[i.ref].push_back( std::make_tuple( i.pos ,i.id )) ;
    }
    // sort by pos 
    for( auto & pair : ref_ranks )
    {
        std::sort(pair.second.begin() , pair.second.end());
    }
    // re-assign rank to contigs
    for(const auto & pair : ref_ranks )
    {
        int rank = 0 ;
        for( const auto & tup : pair.second )
        {
            nodes[std::get<1>(tup)].rank = rank ;
            rank ++ ;
        }
    }
    for( auto & a_edge : edges )
    {
        a_edge.rank = get_rank(a_edge.from,a_edge.to);
    }
}
void prepare_depth() {
    std::set<unsigned int> updated ;
    for( auto & pair : nodes ) {
        pair.second.max_branch_len = 1 ;
        pair.second.max_branch_from_id = pair.second.id;
    }
    for( auto & pair : nodes ) {
        if( pair.second.edges.size() == 1 ) {
            auto up  = pair.second.SendDepth() ;
            for( auto x : up ) updated.insert(x);
        }
    }
    do{
        std::set<unsigned int> new_up ;
        for( auto x : updated ) {
            auto & node = nodes.at(x) ;
            auto up  = node.SendDepth() ;
            for( auto x : up ) new_up.insert(x);
        }
        updated.clear();
        std::swap(new_up,updated);
    }while( ! updated.empty() );
}
void node_class_1()
{
    for(auto & pair : nodes )
    {
        auto & a_node = pair.second ;
        if( a_node.BranchNum() == 1 ) 
            a_node.graph_type = MST_AnalysisNode::GraphType::TipNode ;
        else if ( a_node.BranchNum() == 2 )
            a_node.graph_type = MST_AnalysisNode::GraphType::LinearNode;
        else if ( a_node.BranchNum() < 1 )
            a_node.graph_type = MST_AnalysisNode::GraphType::UnknowNode;
        else 
            a_node.graph_type = MST_AnalysisNode::GraphType::JunctionNode;
    }
}

int oppo_branch_len(const MST_AnalysisNode & nieb , unsigned int id )
{
    int len = 0 ;
    for( const auto & pair : nieb.branch_len_map ){
        if( pair.first != id && pair.second > len ) 
            len = pair.second ;
    }
    return len ;
}

void node_class_2()
{
    for(auto & pair : nodes )
    {
        auto & a_node = pair.second ;
        if( a_node.graph_type != MST_AnalysisNode::GraphType::JunctionNode ) 
            continue ;
        int long_neib = 0 ;
        for( int index : a_node.edges )
        {
            const auto & oppo = get_oppo(index,a_node.id);
            if( oppo_branch_len( oppo , a_node.id ) >= 3 ) long_neib ++ ; 
        }
        if( long_neib >= 3 ) 
            a_node.junction_type = MST_AnalysisNode::JunctionType::LongJunction ;
        else
            a_node.junction_type = MST_AnalysisNode::JunctionType::TipJunction;
    }
}

void node_class_3()
{
    for(auto & pair : nodes )
    {
        auto & a_node = pair.second ;
        if( a_node.graph_type != MST_AnalysisNode::GraphType::JunctionNode ) 
            continue ;
        if( a_node.edges.size() >= 5) 
            a_node.l2_type = MST_AnalysisNode::NodeL2::HighL2 ;
        else
            a_node.l2_type = MST_AnalysisNode::NodeL2::LowL2;
    }
}
void node_class_4()
{
    for(auto & pair : nodes )
    {
        auto & a_node = pair.second ;
        if( a_node.graph_type != MST_AnalysisNode::GraphType::JunctionNode ) 
            continue ;
        int neib_linear = 0 ;
        int neib_tip = 0 ;
        int neib_long = 0 ;
        for( int index : a_node.edges )
        {
            const auto & oppo = get_oppo(index,a_node.id);
            if( oppo.graph_type != MST_AnalysisNode::GraphType::JunctionNode ) 
                neib_linear ++ ;
            else {
                if ( oppo.junction_type == MST_AnalysisNode::JunctionType::TipJunction ) 
                    neib_tip ++ ;
                else 
                    neib_long ++ ;
            }
        }
        if( neib_long == 0 && neib_tip == 0 )
            a_node.nieb_type = MST_AnalysisNode::NeibType::SimpleNieb ;
        else if ( neib_long == 0 && neib_tip != 0 ) 
            a_node.nieb_type = MST_AnalysisNode::NeibType::Short_Long;
        else if ( neib_long != 0 && neib_tip == 0 ) 
            a_node.nieb_type = MST_AnalysisNode::NeibType::Long_Log;
        else 
            a_node.nieb_type = MST_AnalysisNode::NeibType::MaxNieb;
    }
}

void create_output_dir(const std::string & dirname )
{
    if( mkdir ( dirname.c_str() , S_IRUSR|S_IWUSR|S_IXUSR| S_IRGRP | S_IXGRP) != 0  )
    {
        std::cerr<<"FATAL : failed to create directory "<<dirname<<" !\n         exit ... \n";
        exit(1);
    }
    output_dir = dirname ;
}

void report(const std::string &  log)
{
    static std::ostream * report_file = NULL ;
    if( report_file == NULL )
        report_file = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_dir + "/report.txt");
    if( log == "END" && report_file != NULL )
    {
        delete report_file ;
        return ;
    }
    (*report_file) << log <<'\n';
    std::cout<<log<<'\n';
}
bool nied_non_unique(const MST_AnalysisNode & node ) {
    for( int index : node.edges )
        if( ! get_oppo(index,node.id).is_unique )
            return true ;
    return false ;
}

void printNonUnique( std::set<MST_AnalysisNode> & node ) {
    int nn = 0 , nu = 0 , un = 0 , uu = 0 ;
    for( const auto & i : node ) {
        bool nieb_unique =! nied_non_unique(i) ;
        if( i.is_unique && nieb_unique ) uu ++ ;
        if( ! i.is_unique && nieb_unique ) nu ++ ;
        if( i.is_unique && ! nieb_unique ) un ++ ;
        if( ! i.is_unique && !nieb_unique ) nn ++ ;
    }
    report("        # self-non-unique+nieb-non-unique   :" + std::to_string(nn));
    report("        # self-unique+nieb-non-unique       :" + std::to_string(un));
    report("        # self-non-unique+nieb-unique       :" + std::to_string(nu));
    report("        # self-unique+nieb-unique           :" + std::to_string(uu));
}

typedef std::function<bool(const MST_AnalysisNode & node)> NodeFilter ;
void printLeftInfo( const std::string & str ,NodeFilter f )
{
    std::set<MST_AnalysisNode> tmp_s;
    std::set<MST_AnalysisNode> tmp_sl ;
    std::set<MST_AnalysisNode> tmp_ll ;
    std::set<MST_AnalysisNode> tmp_m ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( f(i) ) {
            if( i.nieb_type == MST_AnalysisNode::NeibType::SimpleNieb ) tmp_s.insert(i);
            if( i.nieb_type == MST_AnalysisNode::NeibType::Short_Long) tmp_sl.insert(i);
            if( i.nieb_type == MST_AnalysisNode::NeibType::Long_Log ) tmp_ll.insert(i);
            if( i.nieb_type == MST_AnalysisNode::NeibType::MaxNieb ) tmp_m.insert(i);
        }
    }
    report("    #   "+str +" simple nieb    :" + std::to_string(tmp_s.size()));
    printNonUnique(tmp_s);
    report("    #   "+str +" short_long nieb    :" + std::to_string(tmp_sl.size()));
    printNonUnique(tmp_sl);
    report("    #   "+str +" long long nieb    :" + std::to_string(tmp_ll.size()));
    printNonUnique(tmp_ll);
    report("    #   "+str +" mixed nieb    :" + std::to_string(tmp_m.size()));
    printNonUnique(tmp_m);
}

void print_edge_rank_violin_csv()
{
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_dir + "/edges_rank_violin.csv");
    (*out)<<"id , type , js \n";
    int index = 0 ;
    for( const auto & edge : edges ) 
    {
        (*out)<<++index<<" ,";
        if( edge.rank > 0 &&edge.rank < 4 )
            (*out)<<"rank_"<<edge.rank<<" ,"<<edge.sim<<'\n';
        else if ( edge.rank == 0 )
            (*out)<<"non_unique , "<<edge.sim<<'\n';
        else if ( edge.rank >= 4 )
            (*out)<<"rank_gt_3 ,"<<edge.sim<<'\n';
        else if ( edge.rank == -2 )
            (*out)<<"cross_ref ,"<<edge.sim<<'\n';
        else ;
    }
    delete out;
}
void print_edge_rank_csv()
{
    std::vector<float> rank[6];
    for( const auto & edge : edges ) 
    {
        if( edge.rank < 4  && edge.rank >= 0)
            rank[edge.rank].push_back(edge.sim);
        else if ( edge.rank >= 4 )
            rank[4].push_back(edge.sim);
        else if (edge.rank == -2 )
            rank[5].push_back(edge.sim);
        else 
            assert(0) ;
    }
    for( int i = 0 ; i < 6 ; i++ )
        std::sort(rank[i].rbegin(),rank[i].rend());
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_dir + "/edges_rank.csv");
    (*out)<<"non_unique , rank_1 , rank_2 , rank_3 , rank_gt_3 , rank_cross_ref\n";
    int end_columns = 0 ;
    int index = 0 ;
    while(end_columns < 6 ) 
    {
        end_columns = 0 ;
        for( int i = 0 ; i < 6 ; i ++ )
        {
            if( (int)rank[i].size() > index )
                (*out)<<rank[i][index]<<" , " ;
            else 
            {
                (*out)<<" , " ;
                end_columns ++ ; 
            }
        }
        (*out)<<'\n';
        index ++ ;
    }
    delete out ;
}

/********************************************************
  *
  * check params functions
  *
  ******************************************************/
void check_file_read(const std::string & file_name)
{
    auto t = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file_name);
    if( t == NULL )
    {
        std::cerr<<"FATAL : failed to open "<<file_name<<" to read!"
            <<"\n         exit ... \n";
        exit(1);
    }
    delete t;
}

void check_file_write(const std::string & file_name)
{
    auto t = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file_name);
    if( t == NULL )
    {
        std::cerr<<"FATAL : failed to open "<<file_name<<" to write!"
            <<"\n         exit ... \n";
        exit(1);
    }
    delete t;
}

int subgraphnum()
{
    std::set<unsigned int> reps ;
    for( const auto & pair : nodes )
    {
        reps.insert(dis_joint_set.GetGroup(pair.first));
    }
    return reps.size();
}

/********************************************************
  *
  * main function
  *
  ******************************************************/

int main( int argc , char ** argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,sorted_unique," sorted_unique_contigs.txt file .");
        DEFINE_ARG_REQUIRED(std::string,mintree," mintree.file .");
        DEFINE_ARG_REQUIRED(std::string,output_dir," output directoy name .");
    END_PARSE_ARGS

    check_file_read(sorted_unique.to_string());
    check_file_read(mintree.to_string());
    create_output_dir(output_dir.to_string() );
    report("############    Summary report  #####################\n");
    { // min_tree
        auto mst = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(mintree.to_string());
        report("#  edges                :   "+std::to_string(parse_mintree(*mst)));
        report("#  nodes                :   "+std::to_string(nodes.size()));
        report("#  subgraph             :   "+std::to_string(subgraphnum()));
        delete mst ;
    }
    prepare_depth();
    node_class_1();
    node_class_2();
    node_class_3();
    node_class_4();
    { // sort_unique_contig.txt
        auto sorts = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(sorted_unique.to_string());
        report("# unique-nodes          :   "+std::to_string(parse_sort_unique_contigs(*sorts)));
        delete sorts ;
    }
    auto tjl = [](const MST_AnalysisNode&i){ 
        return i.graph_type == MST_AnalysisNode::GraphType::JunctionNode 
        && i.junction_type == MST_AnalysisNode::JunctionType::TipJunction 
        && i.l2_type == MST_AnalysisNode::NodeL2::LowL2 ;
    }; 
    auto tjh = [](const MST_AnalysisNode&i){ 
        return i.graph_type == MST_AnalysisNode::GraphType::JunctionNode 
        && i.junction_type == MST_AnalysisNode::JunctionType::TipJunction 
        && i.l2_type == MST_AnalysisNode::NodeL2::HighL2 ;
    }; 
    auto ljl = [](const MST_AnalysisNode&i){ 
        return i.graph_type == MST_AnalysisNode::GraphType::JunctionNode 
        && i.junction_type == MST_AnalysisNode::JunctionType::LongJunction
        && i.l2_type == MST_AnalysisNode::NodeL2::LowL2 ;
    }; 
    auto ljh = [](const MST_AnalysisNode&i){ 
        return i.graph_type == MST_AnalysisNode::GraphType::JunctionNode 
        && i.junction_type == MST_AnalysisNode::JunctionType::LongJunction
        && i.l2_type == MST_AnalysisNode::NodeL2::HighL2 ;
    }; 
    prepare_ref_ranks();
        report("# non-unique-nodes      :   "+std::to_string(count_non_unique()));
        report("# junction nodes        :   "+std::to_string(count_junction()));
        report("  ----------------------------------------");
        report("#   junction 3 nodes    :   "+std::to_string(count_junction(3)));
        report("#   junction 4 nodes    :   "+std::to_string(count_junction(4)));
        report("#   junction>4 nodes    :   "+std::to_string(count_junction_gt(4)));
        report("  ----------------------------------------");
        report("#   tip-junction nodes  :   "+std::to_string(count_tip_junction()));
        report("#       low branches    :   "+std::to_string(count_tip_junction_low()));
        printLeftInfo("" ,tjl);
        report("#       high branches   :   "+std::to_string(count_tip_junction_high()));
        printLeftInfo("" ,tjh);
        report("#   long-junction nodes :   "+std::to_string(count_long_junction()));
        report("#       low branches    :   "+std::to_string(count_long_junction_low()));
        printLeftInfo("" ,ljl);
        report("#       high branches   :   "+std::to_string(count_long_junction_high()));
        printLeftInfo("" ,ljh);
        report("  ----------------------------------------");
        report("#   junction unique     :   "+std::to_string(count_junction_unique()));
        report("#   junction non-unique :   "+std::to_string(count_junction_non_unique()));
        report("# tip nodes             :   "+std::to_string(count_tip()));
        report("# linear nodes          :   "+std::to_string(count_linear()));
        report("#   linear nodes 2 rank1 edge :   "+std::to_string(count_linear_edge_rank1(2)));
        report("#   linear nodes 1 rank1 edge :   "+std::to_string(count_linear_edge_rank1(1)));
        report("#   linear nodes 0 rank1 edge :   "+std::to_string(count_linear_edge_rank1(0)));
        report("# rank 1 edge           :   "+std::to_string(count_edge_rank(1)));
        report("# rank 2 edge           :   "+std::to_string(count_edge_rank(2)));
        report("# rank 3 edge           :   "+std::to_string(count_edge_rank(3)));
        report("# rank>3 edge           :   "+std::to_string(count_edge_rank_gt(3)));
        report("# non-unique edge       :   "+std::to_string(count_edge_rank(0)));
        report("# cross-ref edge        :   "+std::to_string(count_edge_rank(-2)));
    report("############    Summary end    #####################\n");
    report("END");
    print_edge_rank_csv();
    print_edge_rank_violin_csv();
    return 0 ;
}
