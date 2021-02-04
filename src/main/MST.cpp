#include "utils/args/argsparser.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"

#include "utils/graph/GraphBasic.h"
#include "utils/graph/mst/MinTree.h"
#include "utils/graph/DisJointSet.h"
#include "utils/graph/GraphTipRemove.h"

#include "stLFR/ContigBinBarcode.h"
#include "utils/misc/fileName.h"

#include <set>
#include <vector>
#include <stack>
#include <map>
#include <sstream>

typedef BGIQD::GRAPH::IGraphNodeBasic<unsigned int , long> Node;

struct Edge : public BGIQD::GRAPH::IGraphEdgeBasic<unsigned int , long >
{
    float sim ;

    std::string AttrString() const
    {
        std::ostringstream ost;
        //ost<<" id= "<<id<<" sim= "<<sim;
        ost<<"label=\""<<sim<<"\"";
        return ost.str();
    }
    std::string ToString() const
    {
        std::ostringstream ost;
        ost<<from<<"\t--\t"<<to<<" [ "<<AttrString()<<" ]";
        return ost.str();
    }
};

// Extract edge weight for MST algorithm
struct EdgeAttr
{
    float GetValue(const Edge & e ) const 
    {
        return 1.0f - e.sim ;
    }
};

struct ContigSimGraph : public BGIQD::GRAPH::Graph<Node,Edge>
{

    typedef BGIQD::GRAPH::Graph<Node,Edge> Basic ;
    typedef Basic::NodeId NodeId;
    typedef Basic::EdgeId EdgeId;


    void AddEdgeSim( unsigned int from , unsigned int to , float sim)
    {
        Edge tmp ;
        tmp.from = from ;
        tmp.to = to ;
        tmp.sim = sim ;
        if( ! HasNode(from) ) {
            Node tnode ;
            tnode.id = from ;
            Basic::AddNode(tnode);
        }
        if( ! HasNode(to) ) {
            Node tnode ;
            tnode.id = to;
            Basic::AddNode(tnode);
        }
        Basic::AddEdge(tmp);
    }

    typedef BGIQD::GRAPH::MinTreeHelper<ContigSimGraph, float , EdgeAttr> MTHelper;
    typedef BGIQD::GRAPH::TipRemoveHelper<ContigSimGraph> TipHelper ;
    typedef TipHelper::TipRemoveResult TipRemoveResult;

    typedef BGIQD::Algorithm::DisJoin_Set<NodeId> DJ_Sets;


    static TipRemoveResult RemoveTip_n2( ContigSimGraph & mintree )
    {
        TipHelper tip_helper ;
        tip_helper.Init(2);
        return tip_helper.DeepTipRemove(mintree);
    }


    static std::vector<NodeId> DetectJunctions( const ContigSimGraph & mintree)
    {
        std::vector<NodeId> ret ;
        for( const auto & pair : mintree.nodes )
        {
            const auto & node = pair.second ;
            if( node.EdgeNum() > 2 )
            {
                ret.push_back( node.id ) ;
            }
        }
        return ret;
    }

    // Get MST by prime-MST algorithm
    ContigSimGraph  MinTree() const 
    {
        EdgeAttr attr ;
        MTHelper helper;
        return helper.MinTree(*this , attr);
    };

    // Get connect componets by disjion set algorithm
    static std::map<NodeId , ContigSimGraph>  UnicomGraph(const ContigSimGraph & mintree)
    {
        std::map<NodeId , ContigSimGraph> ret;
        DJ_Sets dj_sets;
        for( const auto & edge : mintree.edges )
        {
            if( edge.IsValid())
                dj_sets.AddConnect(edge.from , edge.to);
        }
        for( const auto & edge : mintree.edges )
        {
            if( ! edge.IsValid() )
                continue;
            auto rep = dj_sets.GetGroup(edge.from);
            if( ! ret[rep].HasNode(edge.from) )
            {
                ret[rep].AddNode(mintree.GetNode(edge.from));
                ret[rep].GetNode(edge.from).CleanEdges();
            }
            if( ! ret[rep].HasNode(edge.to) )
            {
                ret[rep].AddNode(mintree.GetNode(edge.to));
                ret[rep].GetNode(edge.to).CleanEdges();
            }
            ret[rep].AddEdge(edge);
        }
        return ret ;
    }

    // Print as DOT
    void PrintAsDOT(std::ostream & out) const
    {
        out<<Edge::DOTHead()<<std::endl;
        std::set<int> multis;
        for( const auto & e : edges )
        {
            out<<"\t"<<e.ToString()<<std::endl;
            if( GetNode(e.from).EdgeNum() > 2 )
                multis.insert(e.from);
            if(  GetNode(e.to).EdgeNum()>2 )
                multis.insert(e.to);
        }
        for( int i : multis )
        {
            out<<i<< " [ shape=box]  "<<std::endl;
        }
        out<<"}"<<std::endl;
    }

    // Extract the linear order from a already lineared graph.
    //
    static std::vector<NodeId> TrunkLinear(const ContigSimGraph &  base )
    {
        std::vector<NodeId> ret ;
        NodeId starter ;
        // Find a tip node
        for( const auto & pair : base.nodes )
        {
            auto & node = pair.second ;
            if( node.EdgeNum() == 1 )
            {
                starter = node.id ;
                break;
            }
            else
            {
                assert( node.EdgeNum() == 2 );
            }
        }
        ret.push_back( starter ) ;

        auto & node = base.GetNode(starter);
        ContigSimGraph::Node::NodeEdgeIdIterator begin,end;
        std::tie(begin,end) = node.GetEdges();
        auto edge_id = *begin;
        auto & edge = base.GetEdge(edge_id);
        NodeId next = edge.OppoNode(starter);
        NodeId curr = starter;
        // Travel to the other tip node
        while(1)
        {
            auto & next_node = base.GetNode(next) ;
            ret.push_back(next);
            if( next_node.EdgeNum() == 1 )
            {
                break;
            }
            else
            {
                assert( next_node.EdgeNum() == 2 );
                ContigSimGraph::Node::NodeEdgeIdIterator begin,end;
                std::tie(begin,end) = next_node.GetEdges();
                auto edge_id1 = *(begin);
                auto edge_id2 = *(std::next(begin));
                auto &edge1 = base.GetEdge(edge_id1);
                auto &edge2 = base.GetEdge(edge_id2);
                // Find next step.
                if( edge1.OppoNode(next) != curr)
                {
                    curr = next ;
                    next = edge1.OppoNode(next) ;
                }
                else
                {
                    assert(edge2.OppoNode(next) != curr);
                    curr = next ;
                    next = edge2.OppoNode(next);
                }
            }
        }
        return ret ;
    }
};
//
// Struct to wrap all global variables and functions
//
struct AppConf
{

    typedef std::vector<ContigSimGraph::NodeId> LinearOrder;

    struct MST_correct
    {
        typedef std::vector<ContigSimGraph::NodeId> JunctionResults;

        struct simplify_log
        {

            int base_basic_node_num ;

            int base_node_num ;

            int mst_node_num ;

            ContigSimGraph::TipRemoveResult tip_result;

            JunctionResults junction_result;

            void Init(int basic , int curr )
            {
                base_basic_node_num = basic  ;
                base_node_num = curr ;
                mst_node_num = 0 ;
                tip_result.Init();
                junction_result.clear();
            }
        };

        void Init( const ContigSimGraph & basic , int mr , float df )
        {
            base_contig_sim_graph = basic ;
            max_correct_loop_num = mr ;
            max_del_fac = df ;
            basic_num = base_contig_sim_graph.nodes.size();
        }

        ContigSimGraph base_contig_sim_graph ;

        ContigSimGraph maximum_span_graph ;

        std::stack<simplify_log> logs;

        simplify_log curr_log;

        int basic_num ;

        int max_correct_loop_num ;

        float max_del_fac ;

        void RoundStart()
        {
            curr_log.Init(basic_num , base_contig_sim_graph.nodes.size());
        }

        void RoundEnd()
        {
            logs.push(curr_log);
        }

        // MAKE SURE THIS INTERFACE CALLED ONLY AT THE END OF PROGRAM
        std::vector<JunctionResults> GetAllJunctionResults() const 
        {
            auto log = logs ;
            std::vector<JunctionResults> ret ;
            while( ! log.empty() )
            {
                auto i = log.top();
                ret.push_back(i.junction_result);
                log.pop();
            }
            return ret ;
        }

        void Simplify()
        {
            maximum_span_graph = base_contig_sim_graph.MinTree();
            if ( maximum_span_graph.nodes.size() > 3 )
            {
                curr_log.tip_result = ContigSimGraph
                    ::RemoveTip_n2(maximum_span_graph) ;
                curr_log.junction_result = ContigSimGraph
                    ::DetectJunctions(maximum_span_graph) ;

                if( curr_log.junction_result.size() > 0 )
                {
                    for( auto x : curr_log.junction_result )
                    {
                        base_contig_sim_graph.RemoveNode(x);
                        masked_nodes.insert(x);
                    }
                }
            }
            curr_log.mst_node_num = maximum_span_graph.nodes.size() ;
        }

        bool IsClean() const 
        {
            if ( logs.empty() )
                return false ;

            const auto & item = logs.top() ;

            return item.base_node_num < 4
                || ( int(logs.size()) >= max_correct_loop_num )
                || ( float(item.base_node_num-item.junction_result.size())  / float(item.base_basic_node_num) <= max_del_fac )
                || ( /*item.tip_result.tip_num == 0 &&*/ item.junction_result.size() == 0 ) ;
        }

        void CorrectGraph( )
        {
            do {
                RoundStart();
                Simplify();
                RoundEnd();
                //config.GenerateMinTreeTrunks();
            } while( ! IsClean() ) ;
        }

        std::string log_str() const 
        {
            std::ostringstream ost;
            ost<<curr_log.base_basic_node_num<<'\t';
            ost<<curr_log.base_node_num<<'\t';
            ost<<curr_log.mst_node_num<<'\t';
            ost<<curr_log.junction_result.size();

            return ost.str();
        }

        std::vector<LinearOrder> GenerateLinear()
        {
            std::vector<LinearOrder> ret ;
            if( maximum_span_graph.nodes.size() < 1 )
                return ret ;
            if( curr_log.junction_result.size() > 0 )
            {
                for( auto x : curr_log.junction_result )
                {
                    maximum_span_graph.RemoveNode(x);
                    masked_nodes.insert(x);
                }
            }
            auto splits = ContigSimGraph::UnicomGraph(maximum_span_graph);
            for( auto & pair : splits)
            {
                auto a_line = ContigSimGraph::TrunkLinear(pair.second);
                if( a_line.size() > 1)
                    ret.push_back(a_line);
            }
            return ret ;
        }
        std::set< ContigSimGraph::NodeId > masked_nodes;
    };

    BGIQD::MISC::FileNames fNames ;

    ContigSimGraph graph;

    BGIQD::LOG::logger lger;

    std::map< ContigSimGraph::NodeId , MST_correct>  split_graphs;

    int max_correct_loop_num ;

    float max_del_fac ;

    float smallest ;

    void Init(const std::string & prefix, float f)
    {
        fNames.Init(prefix);
        smallest = f ;
        BGIQD::LOG::logfilter::singleton().get("MST", BGIQD::LOG::loglevel::INFO,lger);
    }

    void LoadContigSimGraph()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.cluster("mst"));
        if( in == NULL )
            FATAL("failed to open xxx.cluster for read!!! ");
        auto parseline = [&]( const std::string & line )
        {
            BGIQD::stLFR::ContigRelation tmp;
            tmp.InitFromString(line);
            for( auto x : tmp.sims )
            {
                if( x.second.simularity >= smallest )
                {
                    graph.AddEdgeSim( tmp.contigId , x.first , x.second.simularity );
                }
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
        lger<<BGIQD::LOG::lstart() << "load contig sim graph done "<<BGIQD::LOG::lend() ;
        lger<<BGIQD::LOG::lstart() << "contig-sim graph nodes : "<<graph.nodes.size()<<BGIQD::LOG::lend() ;
        lger<<BGIQD::LOG::lstart() << "contig-sim graph edges : "<<graph.edges.size()<<BGIQD::LOG::lend() ;
    }

    void SplitGraph()
    {
        auto splits = graph.UnicomGraph(graph);
        for(const auto & pair : splits)
        {
            split_graphs[pair.first].Init(pair.second,max_correct_loop_num,max_del_fac);
        }
        lger<<BGIQD::LOG::lstart() << "split contig sim graph into "<<split_graphs.size()<<" sub graph"<<BGIQD::LOG::lend() ;
    }

    void CorrectGraph()
    {
        for( auto & pair : split_graphs)
        {
            pair.second.CorrectGraph();
            lger<<BGIQD::LOG::lstart() << "CorrectLog\t"<<pair.second.log_str()<<BGIQD::LOG::lend() ;
        }
    }

    void GenerateLinears()
    {
        BGIQD::MISC::Freq<int>  trunk_freq;
        int trunk_count = 1;
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mintreetrunklinear());
        if( out3 == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for write !!! ");
        for( auto & pair : split_graphs)
        {
            auto linears = pair.second.GenerateLinear();
            for(const auto & linear : linears )
            {
                (*out3)<<"---\t"<<trunk_count<<"\t---"<<std::endl;
                trunk_count ++ ;
                trunk_freq.Touch(linear.size());
                for(const auto x : linear )
                {
                    (*out3)<<x<<std::endl;
                }
            }
        }
        lger<<BGIQD::LOG::lstart() << "linear freq is :\n "<<trunk_freq.ToString()<<BGIQD::LOG::lend() ;
        delete out3;
    }

    void PrintJunctionNodes()
    {
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mst_error());
        if( out3 == NULL )
            FATAL(" failed to open xxx.mst_error for write !!! ");
        for(const auto & pair : split_graphs)
        {
            for( auto & i: pair.second.masked_nodes)
            {
                (*out3)<<i<<'\n';
            }
        }
        delete out3;
    }

    void PrintTips()
    {
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mst_tips());
        if( out3 == NULL )
            FATAL(" failed to open xxx.mst_tips for write !!! ");
        for(const auto & pair : split_graphs)
        {
            for( auto & i: pair.second.curr_log.tip_result.tip_contigs)
            {
                (*out3)<<i<<'\n';
            }
        }
        delete out3;
    }
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.mintree_trunk_linear && xxx.mst_error && xxx.mst_tips ");
        DEFINE_ARG_OPTIONAL(float , threshold, " threshold of simularity","0.1");
        DEFINE_ARG_OPTIONAL(float , del_fac, " threshold of del junction node","0.9");
        DEFINE_ARG_OPTIONAL(int, del_round, "maximum del round ","1000");
    END_PARSE_ARGS

    config.Init(prefix.to_string(),threshold.to_float());
    config.max_del_fac = del_fac.to_float();
    config.max_correct_loop_num = del_round.to_int() ;
    config.LoadContigSimGraph();
    config.SplitGraph();
    config.CorrectGraph();
    config.GenerateLinears();
    config.PrintJunctionNodes() ;
    config.PrintTips();
    return 0 ;
}
