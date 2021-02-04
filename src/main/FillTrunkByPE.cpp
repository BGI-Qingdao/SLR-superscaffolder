/**********************************************************
 *
 *
 *
 *
 *
 * *******************************************************/
#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/string/stringtools.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"
#include "utils/misc/fileName.h"
#include "stLFR/contigPEGraph.h"

#include <map>
#include <set>
#include <stack>
#include <sstream>
#include <string>
#include <iostream>
#include <cassert>
#include <functional>

using namespace BGIQD::stLFR;


/*********************************************************
 *
 * We wrote an customized shortest-path-search algorithm 
 * based on depth-first search to support an special feature.
 *
 * Awkwardly, the feature that need to design this customized
 * algorithm was abandoned now.
 *
 * Therefor, we will replace this algorithm with a standard
 * shortest-path-search algorithm like Dijkstra.
 *
 * ******************************************************/
enum DepthSearchEdgeType
{
    Invalid = -1 ,
    White = 0 ,
    Gray = 1 ,
    Black = 2 ,
    EndPoint  = 3 ,
};

// The node for searched path graph, follow links of prev can travel to path.
template<class Node>
struct DepthSearchNode
{
    typedef typename Node::NodeNodeId   NodeId;
    typedef DepthSearchEdgeType         Type;
    Type                                type ;

    NodeId                              id;
    NodeId                              prev ;

    int                                 first_found ;
    int                                 finish_search ;

    static const  NodeId                invalid = -1;

    std::set<NodeId>                    backword_from ;
    std::set<NodeId>                    forword_from ;
    std::set<NodeId>                    cross_from ;

    void Clean()
    {
        type = White ;
        prev = invalid ;
        first_found = -1 ;
        finish_search = -1 ;

        backword_from.clear();
        forword_from.clear();
        cross_from.clear();
    }

    void ReSetParent(const Node & me,const DepthSearchNode & parenet,int step_start )
    {
        prev = parenet.id ;
        first_found = step_start ;
        type = White ;
        id = me.id ;
        //TODO : reset edges
    }

    void InitRoot(const Node & me,int step_start )
    {
        type = White ;
        prev = invalid;
        id = me.id;
        first_found =  step_start ;
    }
    void Init(const Node & me,const DepthSearchNode & parenet,int step_start )
    {
        type = White ;
        prev = parenet.id;
        id = me.id;
        first_found =  step_start ;
    }

    std::string ToString() const {
        std::ostringstream ost;
        ost<<"id = "<<id<<" prev= "<<prev;
        ost<<" first= "<<first_found<<" last= "<<finish_search;
        for( auto x : backword_from)
        {
            ost<<" Backword= "<<x;
        }
        for( auto x : forword_from)
        {
            ost<<" Forword= "<<x;
        }
        for( auto x : cross_from)
        {
            ost<<" Cross= "<<x;
        }
        return ost.str();
    }
};

// Design an special edge iterator because the edge should be iterated after sorted by PE linkage weight.
template<class Graph>
struct NodeEdgeIterator
{
    public:

        typedef typename Graph::EdgeId        Id;
        typedef typename Graph::Node          Node;
        typedef typename Graph::Edge          Edge;

        NodeEdgeIterator() : curr(NULL) ,accessor( NULL ) { }

        NodeEdgeIterator(const Node & n,const Edge & e , Graph & acc )
        {
            node_id = n.id ;
            if( e.id != Edge::invalid )
                curr = &e;
            else
                curr = NULL ;
            accessor = &acc ;
        }

        NodeEdgeIterator( const NodeEdgeIterator & ei )
        {
            node_id = ei.node_id;
            curr = ei.curr ;
            accessor = ei.accessor ;
        }

        NodeEdgeIterator & operator = ( const NodeEdgeIterator & ei )
        {
            if( &ei != this )
            {
                node_id = ei.node_id ;
                curr = ei.curr ;
                accessor = ei.accessor ;
            }
            return *this;
        }

        // ONLY detect curr.
        // ONLY use == with End() .
        bool operator == ( const NodeEdgeIterator & ei )const
        {
            return curr == ei.curr ;
        }

        // ONLY detect curr.
        // ONLY use == with End() .
        bool operator != ( const NodeEdgeIterator & ei )const
        {
            return curr != ei.curr ;
        }

        const Edge & operator*() const  { return *curr ; }

        const Edge * operator->() const  { return curr ; }

        NodeEdgeIterator & operator ++() {
            if( curr != NULL && accessor != NULL )
            {
                Id next = curr->next ;
                if( next != Edge::invalid )
                    curr = &(accessor->GetEdge(next));
                else
                    curr = NULL ;
            }
            else
                curr = NULL ;
            return *this ;
        }

        static NodeEdgeIterator & end() 
        {
            static NodeEdgeIterator end;
            return end ;
        }
        typename Edge::EdgeNodeId node_id ;
    private:
        const Edge * curr ;
        Graph      * accessor ;
};
template<class Graph /* the graph to be searched*/
, class EdgeItr      /* an iterator to visit all edge for a node*/
, class PathEnder    /* detecter to give up a search*/
, class DepthNode = DepthSearchNode<typename Graph::Node > >
struct DepthSearch
{
    typedef typename Graph::NodeId   NodeId;
    typedef typename Graph::EdgeId   EdgeId;
    typedef typename Graph::Node          NodeBase;
    typedef typename Graph::Edge          EdgeBase;
    typedef DepthNode                           Node;
    //              typedef DepthSearchEdge<EdgeBase>           Edge;

    NodeId                                      start;

    std::map<NodeId , Node>                     nodes;
    //std::map<EdgeId , Edge>                   edges;
    std::stack<EdgeItr>                         path;
    Graph                                       *accesser;
    PathEnder                                   ender;

    void PrintNodes() const
    {
        for( auto & i : nodes )
        {
            std::cerr<<i.second.ToString()<<std::endl;
        }
    }

    //  do depth search.
    //      1. stop go deeper if ender say yes.
    //      2. avoid to use recursion because the depth may very large.
    //
    // @params
    //          start : id of root node
    //          step  : index of start step
    // @return 
    //          end step index
    //
    int DoDepthSearch(NodeId start , int s_step )
    {
        int step = s_step;
        NodeBase & root = accesser->GetNode(start);
        Node & curr = nodes[start];

        curr.InitRoot(root,step);

        path.push(EdgeItr(root,accesser->GetEdge(root.edge_id ) , *accesser));
        ender.Start();
        bool new_node_in_path = true ;
        NodeId prev ;
        while ( ! path.empty() )
        {
            auto & itr = path.top() ;
            Node * test ;
            try {
                Node & top = nodes.at( itr.node_id ) ;
                test = &top;
            }
            catch(...)
            {
                assert(0);
            }
            Node & top = *test;
            if( itr == EdgeItr::end() )
            {
                step ++ ;
                if( top.type == Node::Type::White )
                    top.type = Node::Type::EndPoint ;
                else
                    top.type = Node::Type::Black ;
                top.finish_search = step ;
                path.pop() ;
                ender.PopEdge();
                if( ! new_node_in_path )
                    ender.PopNode();

                new_node_in_path = false ;
                //assert(path.size() == ender.nodes.size() );
                continue ;
            }

            if( new_node_in_path )
            {
                if( ender.IsEnd() )
                { // end by edge
                    top.Clean();
                    ender.PopEdge();
                    path.pop() ;
                    new_node_in_path = false ;
                    //assert(path.size() == ender.nodes.size() );
                    continue ;
                }

                ender.AddNode(accesser->GetNode(top.id),top);
                if( ender.IsEnd() )
                {
                    if( top.type == Node::Type::White )
                    { // end by new node
                        top.type = Node::Type::EndPoint ;
                        step ++ ;
                        top.finish_search = step ;
                    }
                    else 
                    { // end by old node
                        ;
                    }
                    ender.PopEdge();
                    ender.PopNode();
                    path.pop() ;
                    new_node_in_path = false ;
                    //assert(path.size() == ender.nodes.size() );
                    continue ;
                }

                if( top.type != Node::Type::White )
                {
                    step ++ ;
                    try {
                        top.ReSetParent( accesser->GetNode(top.id) , nodes.at(prev), step );
                    }
                    catch(...)
                    {
                        assert(0);
                    }
                }
                //in path now , force make it Gray !!!
                top.type = Node::Type::Gray ;
            }

            //assert(path.size() == ender.nodes.size() );
            NodeId next_node_id  = itr->to;
            auto itr_n = nodes.find( next_node_id ) ;
            if ( itr_n == nodes.end() || itr_n->second.type == Node::Type::White )
            {
                step ++ ;
                // white node alloc 
                auto & next_node = nodes[next_node_id] ;
                next_node.Init( accesser->GetNode(next_node_id) , top, step );
            }
            else
            {
                auto & node = itr_n->second ;
                if ( node.type == Node::Type::Gray )
                { // match a backwoard edge
                    node.backword_from.insert( top.id ) ;
                }
                else
                {// black node or endpoint left

                    if( node.first_found > top.first_found  )
                    { // to child ,match a forward edge
                        node.forword_from.insert(top.id);
                    }
                    else
                    { // to other branch ,match a corss edge
                        node.cross_from.insert(top.id);
                    }
                }
            }
            prev = top.id ;
            auto & next_node = accesser->GetNode( next_node_id) ;
            EdgeId next_edge_id = next_node.edge_id;
            EdgeBase & next_edge = accesser->GetEdge( next_edge_id) ;
            //assert(path.size() == ender.nodes.size() );
            path.push(EdgeItr(next_node,next_edge, *accesser));
            ender.AddEdge(next_edge);

            new_node_in_path = true ;
            ++ itr ;
        }
        return step ;
    }
};

//
// Give up a search when
//  1. path too depth or too long
//  2. searched a already in ordered-chain contig.
//
struct DepthEnder
{
    bool ender_flag ;
    typedef ContigPEGraph::NodeId NodeId;
    typedef ContigPEGraph::Node   Node;
    typedef ContigPEGraph::Edge   Edge;
    typedef DepthSearchNode<Node> SNode;

    enum NodeType
    {
        Unknow = 0 ,
        Others = 1 ,
        Seed = 2 ,
    };

    typedef std::function<NodeType(const NodeId &)> NodeTypeDetector;

    NodeTypeDetector keyer ;

    int max_length ;

    int curr_len ;

    int last_len ;

    bool first_in ;

    std::stack<const Node * > path;

    void AddNode(const Node & node , const SNode & snode )
    {
        path.push(&node);

        if( ! first_in )
        {
            first_in = true ;
            return ;
        }
        if( snode.type != DepthSearchEdgeType::White )
        {
            ender_flag = true ;
            return ;
        }
        auto ret =  keyer(node.id);
        if( ret == NodeType::Unknow )
        {
            assert(0);
        }
        if( ret == NodeType::Others) 
        {
            if ( max_length != -1 && curr_len > max_length )
            {
                ender_flag = true ;
                return ;
            }
            else
            {
                curr_len += node.contigLen;
                return ;
            }
        }
        else if ( ret == NodeType::Seed)
        {
            ender_flag = true ;
            return ;
        }
        else
            assert(0);
        ender_flag = true ;
    }

    void PopEdge() { ender_flag  = false ; }

    void PopNode()
    {
        if( path.empty() )
        {
            assert(0);
            return ;
        }
        const Node * top = path.top();
        path.pop();
        ender_flag  = false ;
        curr_len = top->contigLen ;
    }

    void AddEdge(const Edge & ) { ender_flag = false ; }

    void Start()
    {
        ender_flag = false ;
        curr_len =last_len = 0 ;
    }
    void Init( NodeTypeDetector  k , int max_l )
    {
        keyer= k ;
        max_length = max_l;
        first_in = false ;
    }
    bool IsEnd() const { return ender_flag ; }
};

//
// Struct to wrap all global variables and functions
//
struct AppConfig
{
    struct GapInfo
    {
        unsigned int prev ;
        unsigned int next ;
        unsigned int prev_1 ;
        unsigned int next_1 ;
        std::set<unsigned int> relations;
        void InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            ist>>prev>>next;
            relations.insert(prev);
            relations.insert(prev+1);
            relations.insert(next);
            relations.insert(next+1);
            while( !ist.eof() )
            {
                unsigned int r ;
                ist>>r ;
                relations.insert(r);
                relations.insert(r+1);
            }
        }
        std::vector<unsigned int> path;
    };

    BGIQD::LOG::logger loger;
    BGIQD::MISC::FileNames fName;
    int searchMax ;
    int total_fill ; 
    std::set<unsigned int> trunk_seeds;
    std::vector<GapInfo> all_gaps;
    BGIQD::stLFR::ContigPEGraph pe_graph;
    std::map<unsigned int , int > contigLen_cache;
    BGIQD::MISC::Freq<int> fill_freq;
    int is_max;
    int min_count ;

    void LoadSeedCluster()
    {
        auto parseline = [this](const std::string & line) -> void 
        {
            GapInfo tmp;
            tmp.InitFromString(line);
            trunk_seeds.insert(tmp.prev);
            trunk_seeds.insert(tmp.next);
            all_gaps.push_back(tmp);
        };
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds_cluster_seeds());
        if( in == NULL )
            FATAL( " failed to open xxx.seeds_cluster_seeds to read !!!" );

        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
    }

    void LoadPEGraph()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.pe_graph());
        if( in == NULL )
            FATAL( " failed to open xxx.pe_graph to read !!!" );

        auto parseline = [this](const std::string & line) -> void 
        {
            if( line == "digraph {" || line == "}" )
                return ;
            BGIQD::stLFR::PEEdge edge;
            edge.InitFromString(line);
            if( std::abs(edge.len) > is_max )
                return ;
            if( edge.count < min_count ) 
                return ;
            pe_graph.AddEdge(edge);
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
    }

    void CalcAll()
    {

        auto typer = [this](unsigned int id )
        {
            if( trunk_seeds.find(id) != trunk_seeds.end() )
                return DepthEnder::NodeType::Seed ;
            else
                return DepthEnder::NodeType::Others;
        };

        typedef NodeEdgeIterator<BGIQD::stLFR::ContigPEGraph> EdgeItr;

        typedef DepthSearch<
            BGIQD::stLFR::ContigPEGraph,
            EdgeItr,
            DepthEnder
                > Searcher;

        auto get_path = [](const Searcher & searcher , unsigned int root , unsigned int end )
        {
            std::vector<unsigned int> ret ;
            if( searcher.nodes.find( end) == searcher.nodes.end() )
                return ret;
            std::stack<unsigned int> path;
            auto & ender = searcher.nodes.at( end ) ;
            unsigned int prev = ender.prev ;
            int max_depth = 1000 ;
            int curr = 0 ;
            path.push( end );
            while( prev != root  && curr < max_depth )
            {
                path.push(prev) ;
                auto & prever = searcher.nodes.at(prev) ;
                prev = prever.prev ;
                curr ++ ;
            }
            if( prev== root )
            {
                path.push(root);
                while(!path.empty() )
                {
                    ret.push_back( path.top() );
                    path.pop();
                }
            }
            return ret;
        };

        for(auto & node : pe_graph.nodes )
        {
            pe_graph.MakeEdgeNext(node.second.id);
        }

        for (auto & gap : all_gaps)
        {
            BGIQD::stLFR::ContigPEGraph cluster_sub_graph;
            BGIQD::stLFR::ContigPEGraph * sub_graph_ptr ;
            cluster_sub_graph = pe_graph.SubGraph(gap.relations);
            sub_graph_ptr = &cluster_sub_graph ;
            BGIQD::stLFR::ContigPEGraph & sub_graph = *(sub_graph_ptr);
            for(auto & node : sub_graph.nodes )
            {
                sub_graph.MakeEdgeNext(node.second.id);
            }
            std::vector<std::vector<unsigned int> > paths;
            unsigned int s1 , e1 ;
            //try 1 order
            {
                Searcher searcher;
                searcher.accesser = &sub_graph;
                searcher.ender.Init(typer, searchMax );
                searcher.DoDepthSearch(gap.prev,0);
                if( searcher.nodes.find(gap.next) != searcher.nodes.end()  )
                {
                    s1= gap.prev ;
                    e1 = gap.next;
                    paths.push_back(get_path(searcher,gap.prev,gap.next));
                }
                if( searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    paths.push_back(get_path(searcher,gap.prev,gap.next+1));
                    s1= gap.prev ;
                    e1 = gap.next +1;
                }
            }
            // try another order
            {
                typedef NodeEdgeIterator<BGIQD::stLFR::ContigPEGraph> EdgeItr;

                typedef DepthSearch<
                    BGIQD::stLFR::ContigPEGraph,
                    EdgeItr,
                    DepthEnder
                        > Searcher;
                Searcher searcher;
                searcher.accesser = &sub_graph;
                searcher.ender.Init(typer, searchMax );
                searcher.DoDepthSearch(gap.prev+1,0);
                if( searcher.nodes.find(gap.next) != searcher.nodes.end()  )
                {
                    paths.push_back(get_path(searcher,gap.prev+1,gap.next));
                    s1= gap.prev +1 ;
                    e1 = gap.next;
                }
                if( searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    paths.push_back(get_path(searcher,gap.prev+1,gap.next+1));
                    s1= gap.prev + 1 ;
                    e1 = gap.next +1;
                }
            }
            fill_freq.Touch(paths.size());
            if( paths.size() == 1 )
            {
                gap.path = paths[0] ;
                gap.prev_1 = s1 ;
                gap.next_1 = e1 ;
            }
        }
        loger<<BGIQD::LOG::lstart()<<fill_freq.ToString()<<BGIQD::LOG::lend();
    }

    void PrintFillResult()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.trunk_fill());
        if( out == NULL )
            FATAL(" failed to open xxx.trunk_fill to write !! ");
        for( const auto & gap :all_gaps )
        {
            if( gap.path.empty() )
                continue ;

            (*out)<<gap.prev<<'\t'<<gap.next<<'\t'
                <<gap.prev_1<<'\t'<<gap.next_1;
            for( auto i : gap.path )
            {
                (*out)<<'\t'<<i;
            }
            (*out)<<'\n';
        }
        delete out ;
    }

    void Init(const std::string & prefix , int search_len_max , int max_is )
    {
        fName.Init(prefix);
        searchMax = search_len_max ;
        BGIQD::LOG::logfilter::singleton().get("FillTrunkByPE", BGIQD::LOG::loglevel::INFO,loger);
        total_fill = 0 ;
        is_max = max_is ;
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds("pe")) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            auto items = BGIQD::STRING::split(line,"\t");
            contigLen_cache[std::stoul(items[0])] = std::stoul(items[1]);
            pe_graph.AddNode(std::stoul(items[0]) , std::stoul(items[1]));
        }
        delete in ;
    };

} config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files.");
        DEFINE_ARG_OPTIONAL(int, searchMax,"max search length." ,"5000");
        DEFINE_ARG_OPTIONAL(int, insert_max,"max insert_size avaliable." ,"1000");
        DEFINE_ARG_OPTIONAL(int, min_count ,"min valid count ." ,"2");
    END_PARSE_ARGS;
    config.Init(prefix.to_string() , searchMax.to_int(),insert_max.to_int());
    config.min_count = min_count.to_int();
    BGIQD::LOG::timer t(config.loger,"FillTrunkByPE");
    config.LoadSeedCluster();
    config.LoadSeeds();
    config.LoadPEGraph();
    config.CalcAll();
    config.PrintFillResult();
    return 0;
}
