/**********************************************************
  *
  *
  *
  *
  *
  ********************************************************/

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
#include <vector>
#include <tuple>

/********************************************************
  *
  * basic structures
  *
  ******************************************************/

struct MST_AnalysisEdge {
    unsigned int from ;
    unsigned int to ;
    int rank ; // 0 means non-unique relative edge
    float sim ;
};

struct MST_AnalysisNode {

    unsigned int id ;

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
        BrankJunction = 2 ,
        MixedJunction = 3 ,
        NeibJunction = 4 
    };

    JunctionType junction_type ;

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
        sscanf(line.c_str(), "%llu\t--\t%llu [ lable=\"%f\" ]",&from, &to,&sim);
        edges.push_back({ (unsigned int )from , (unsigned int)to , 0 , sim});
        nodes[from].AddEdge(from,edges.size()-1);
        nodes[to].AddEdge(to ,edges.size()-1);
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

const MST_AnalysisNode & get_oppo( int edge_index , unsigned int self_id )
{
    const auto & the_edge = edges.at(edge_index) ;
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
    if( from_node.ref != to_node.ref )  return 0 ;
    return std::abs(from_node.rank  - to_node.rank );
}

int count_edge_rank(int num)
{
    std::string file ;
    if( num != 0 )
        file = output_dir+"/edge_rank_"+std::to_string(num)+".txt";
    else
        file = output_dir+"/edge_non_unique.txt";
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
        if( i.junction_type == MST_AnalysisNode::JunctionType::BrankJunction  ) ret ++ ;
    }
    return ret ;
}
int count_long_junction_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( ! i.is_unique && i.junction_type == MST_AnalysisNode::JunctionType::BrankJunction  ) ret ++ ;
    }
    return ret ;
}
int count_long_junction_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::BrankJunction  ) 
        {
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}
int count_long_junction_or_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::BrankJunction  ) 
        {
            if( ! i.is_unique )
            {
                ret ++ ;
                continue ;
            }
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}


int count_mixed_junction()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::MixedJunction) ret ++ ;
    }
    return ret ;
}
int count_mixed_junction_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( ! i.is_unique && i.junction_type == MST_AnalysisNode::JunctionType::MixedJunction) ret ++ ;
    }
    return ret ;
}
int count_mixed_junction_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::MixedJunction) 
        {
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}
int count_mixed_junction_or_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::MixedJunction) 
        {
            if( ! i.is_unique )
            {
                ret ++ ;
                continue ;
            }
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
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
int count_tip_junction_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( ! i.is_unique && i.junction_type == MST_AnalysisNode::JunctionType::TipJunction  ) ret ++ ;
    }
    return ret ;
}
int count_tip_junction_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::TipJunction  ) 
        {
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}
int count_tip_junction_or_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::TipJunction  ) 
        {
            if( ! i.is_unique )
            {
                ret ++ ;
                continue ;
            }
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}

int count_nieb_junction() 
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( i.junction_type == MST_AnalysisNode::JunctionType::NeibJunction) ret ++ ;
    }
    return ret ;
}
int count_nieb_junction_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if( ! i.is_unique && i.junction_type == MST_AnalysisNode::JunctionType::NeibJunction) ret ++ ;
    }
    return ret ;
}
int count_nieb_junction_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::NeibJunction) 
        {
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
    }
    return ret ;
}
int count_nieb_junction_or_nib_non_unique()
{
    int ret = 0 ;
    for( const auto & pair : nodes)
    {
        const auto & i = pair.second ;
        if(  i.junction_type == MST_AnalysisNode::JunctionType::NeibJunction) 
        {
            if( ! i.is_unique )
            {
                ret ++ ;
                continue ;
            }
            for( int index : i.edges )
                if( ! get_oppo(index,i.id).is_unique )
                {
                    ret ++ ;
                    break ;
                }
        }
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


struct JunctionTypeDetecter {
    int tip ;
    int linear ;
    int junction ;
    JunctionTypeDetecter() : tip(0) , linear(0) , junction(0) {}
    void Add( MST_AnalysisNode::GraphType t)
    {
        if( t == MST_AnalysisNode::GraphType::JunctionNode ) junction ++ ;
        else if ( t == MST_AnalysisNode::GraphType::LinearNode ) linear ++ ;
        else if ( t == MST_AnalysisNode::GraphType::TipNode ) tip ++ ;
        else assert(0);
    }
    MST_AnalysisNode::JunctionType get_type() const 
    {
        if( junction > 0 ) return MST_AnalysisNode::JunctionType::NeibJunction ;
        if( linear <= 2 && tip >= 1 && junction == 0 ) return MST_AnalysisNode::JunctionType::TipJunction ;
        if( linear > 2 && tip == 0 && junction == 0 ) return MST_AnalysisNode::JunctionType::BrankJunction ;
        if( linear > 2 && tip > 0 && junction == 0 ) return MST_AnalysisNode::JunctionType::MixedJunction ;
        assert(0);
        return MST_AnalysisNode::JunctionType::UnknowJunction ;
    }
};
void node_class_2()
{
    for(auto & pair : nodes )
    {
        auto & a_node = pair.second ;
        if( a_node.graph_type != MST_AnalysisNode::GraphType::JunctionNode ) 
            continue ;
        JunctionTypeDetecter detecter;
        for( int index : a_node.edges )
        {
            const auto & oppo = get_oppo(index,a_node.id);
            detecter.Add(oppo.graph_type);
        }
        a_node.junction_type = detecter.get_type();
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

void print_edge_rank_csv()
{
    std::vector<float> rank[5];
    for( const auto & edge : edges ) 
    {
        if( edge.rank < 4 )
            rank[edge.rank].push_back(edge.sim);
        else
            rank[4].push_back(edge.sim);
    }
    for( int i = 0 ; i < 5 ; i++ )
        std::sort(rank[i].rbegin(),rank[i].rend());
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_dir + "/edges_rank.csv");
    (*out)<<"non-unique , rank-1 , rank-2 , rank-3 , rank>4\n";
    int end_columns = 0 ;
    int index = 0 ;
    while(end_columns < 5 ) 
    {
        end_columns = 0 ;
        for( int i = 0 ; i < 5 ; i ++ )
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
        delete mst ;
    }
    node_class_1();
    node_class_2();
    { // sort_unique_contig.txt
        auto sorts = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(sorted_unique.to_string());
        report("# unique-nodes          :   "+std::to_string(parse_sort_unique_contigs(*sorts)));
        delete sorts ;
    }
    prepare_ref_ranks();
        report("# non-unique-nodes      :   "+std::to_string(count_non_unique()));
        report("# junction nodes        :   "+std::to_string(count_junction()));
        report("  ----------------------------------------");
        report("#   junction 3 nodes    :   "+std::to_string(count_junction(3)));
        report("#   junction 4 nodes    :   "+std::to_string(count_junction(4)));
        report("#   junction>4 nodes    :   "+std::to_string(count_junction_gt(4)));
        report("  ----------------------------------------");
        report("#   tip-junction nodes  :   "+std::to_string(count_tip_junction()));
        report("#     tip-junction nodes or nibs non-unique :   "+std::to_string(count_tip_junction_or_nib_non_unique()));
        report("#     tip-junction nodes non-unique         :   "+std::to_string(count_tip_junction_non_unique()));
        report("#     tip-junction nodes nibs non-unique    :   "+std::to_string(count_tip_junction_nib_non_unique()));
        report("      ------------------------------------");
        report("#   long-junction nodes :   "+std::to_string(count_long_junction()));
        report("#     long-junction nodes or nibs non-unique:   "+std::to_string(count_long_junction_or_nib_non_unique()));
        report("#     long-junction nodes non-unique        :   "+std::to_string(count_long_junction_non_unique()));
        report("#     long-junction nodes nibs non-unique   :   "+std::to_string(count_long_junction_nib_non_unique()));
        report("      ------------------------------------");
        report("#   mixed-junction nodes:   "+std::to_string(count_mixed_junction()));
        report("#     mixed-junction nodes or nibs non-unique:   "+std::to_string(count_mixed_junction_or_nib_non_unique()));
        report("#     mixed-junction nodes non-unique        :   "+std::to_string(count_mixed_junction_non_unique()));
        report("#     mixed-junction nodes nibs non-unique   :   "+std::to_string(count_mixed_junction_nib_non_unique()));
        report("      ------------------------------------");
        report("#   nieb-junction nodes :   "+std::to_string(count_nieb_junction()));
        report("#     nieb-junction nodes or nibs non-unique:   "+std::to_string(count_nieb_junction_or_nib_non_unique()));
        report("#     nieb-junction nodes non-unique        :   "+std::to_string(count_nieb_junction_non_unique()));
        report("#     nieb-junction nodes nibs non-unique   :   "+std::to_string(count_nieb_junction_nib_non_unique()));
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
    report("############    Summary end    #####################\n");
    report("END");
    print_edge_rank_csv();
    return 0 ;
}
