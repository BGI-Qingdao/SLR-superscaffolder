#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "stLFR/CBB.h"
#include "stLFR/contigSimGraph.h"
#include "stLFR/CBB.h"

#include "soap2/fileName.h"

#include <set>
#include <vector>
#include <stack>
#include <tuple>
#include <sstream>

struct AppConf
{

    typedef std::vector<BGIQD::stLFR::ContigSimGraph::NodeId> LinearOrder;

    struct MST_correct_new 
    {
        struct G1Edge : public BGIQD::GRAPH::IGraphEdgeBasic<unsigned int , int >
        {
            int weight;
            void IncrWeight() { weight ++ ; } 
        };

        struct G1Node : public BGIQD::GRAPH::IGraphNodeBasic<unsigned , int >
        {

        };

        struct GraphG1 : public BGIQD::GRAPH::ListGraph<G1Node,G1Edge> 
        {
            typedef BGIQD::Algorithm::DisJoin_Set<NodeId> DJ_Sets;
            std::map<unsigned int , GraphG1 > UnicomGraph()  const 
            {
                std::map<NodeId , GraphG1 > ret;
                DJ_Sets dj_sets;
                for( const auto & edge : edges )
                {
                    if( edge.IsValid())
                        dj_sets.AddConnect(edge.from , edge.to);
                }
                for( const auto & edge : edges )
                {
                    if( ! edge.IsValid() )
                        continue;
                    auto rep = dj_sets.GetGroup(edge.from);
                    if( ! ret[rep].HasNode(edge.from) )
                    {
                        ret[rep].AddNode(GetNode(edge.from));
                        ret[rep].GetNode(edge.from).edge_ids.clear();
                    }
                    if( ! ret[rep].HasNode(edge.to) )
                    {
                        ret[rep].AddNode(GetNode(edge.to));
                        ret[rep].GetNode(edge.to).edge_ids.clear();
                    }
                    ret[rep].AddEdge(edge);
                }
                return ret ;
            }

            typedef BGIQD::GRAPH::ListGraph<G1Node,G1Edge> Basic;
            void InitEdge( unsigned int from , unsigned int to )
            {
                Edge tmp ;
                tmp.from = from ;
                tmp.to = to ;
                tmp.weight = 1;
                Basic::AddEdge(tmp);
            }
            void Add1Edge(unsigned int from , unsigned int to )
            {
                if( Basic::CheckEdge(from, to) ) {
                    auto & edge = GetEdge(from,to); 
                    edge.IncrWeight();
                } else {
                    InitEdge(from,to);
                }
            }
        };

        struct OldRoad {
            unsigned int left ;
            unsigned int mid ;
            unsigned int right ;
        } ;

        struct NewRoad : public OldRoad {
            bool valid ;

        };

        std::vector<OldRoad> GetCandidateRoads(  const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info ) {
            std::vector<OldRoad> ret ;
            for( int i = 0 ; i < (int)junction_info.neibs.size() ; i ++ )
            {
                for( int j = 0 ; j < (int)junction_info.neibs.size() ; j ++ )
                    ret.emplace_back( OldRoad{ junction_info.neibs.at(i) 
                            , junction_info.junction_id
                            , junction_info.neibs.at(j) } );

            }
            return ret ;
        }

        NewRoad CheckRoad( const OldRoad & old_road ){
            NewRoad ret ;
            ret.valid = false ;
            //TODO
            return ret ;
        }

        bool ContainCircle( const GraphG1 & graph){
            bool ret = false ;
            return ret ;
        }

        GraphG1 GetG1( const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info )
        {
            GraphG1 ret ;
            auto cands = GetCandidateRoads(junction_info);
            for( const auto & cand : cands )
            {
                auto new_road = CheckRoad(cand);
                if( new_road.valid ) {
                    ret.Add1Edge(new_road.left, new_road.mid ) ;
                    ret.Add1Edge(new_road.right, new_road.mid ) ;
                }
            }
            if( ! ContainCircle(ret ) ) {
            for( auto neib : junction_info.neibs )
                ret.Add1Edge(neib,junction_info.junction_id);
            }
            return ret ;
        }

        BGIQD::stLFR::ContigSimGraph base_contig_sim_graph ;
        BGIQD::stLFR::ContigSimGraph mst_v1;
        BGIQD::stLFR::ContigSimGraph mst_v2;
        int min_common_barcode_type ;
        float min_53 ; 

        // for debug print 
        int base_num ; 
        int mst_num ;
        int remain_junction_num ;
        int remain_tip_num ;
        int in_linear_num ;


        void Init(const BGIQD::stLFR::ContigSimGraph & base , int mc , float m53 ) {
            base_contig_sim_graph = base ;
            min_53 = m53 ;
            min_common_barcode_type = mc ;
            base_num = base.NodesSize() ;
        }
        static std::pair< int , GraphG1 > GetA_group( const GraphG1 & graph_g1 ) {
            auto splits = graph_g1.UnicomGraph() ;
            GraphG1 * the_one = NULL ;
            int group_num = 0 ;
            if( splits.size() > 1 ) {
                for( auto & pair : splits ) {
                    if( pair.second.NodesSize() > 1 ) {
                        group_num ++ ;
                        the_one = &pair.second ;
                    }
                }
            }
            assert(group_num > 0 && the_one != NULL );
            return std::make_pair( group_num , *the_one ) ;
        }

        bool CanSimplify( const  GraphG1 & graph_g1 ) {
            auto gret = GetA_group(graph_g1) ;
            if( gret.first > 1  ) return false ;
            const auto & tmp_g1 = gret.second ;
            for( const auto & pair : tmp_g1.nodes ) {
                if( pair.second.EdgeNum() > 2 )
                    return true ;
            }
            return false;
        }

        bool IsLinear( const  GraphG1 & graph_g1) {
            auto gret = GetA_group(graph_g1) ;
            if( gret.first > 1  ) return false ;
            const auto & tmp_g1 = gret.second ;
            for( const auto & pair : tmp_g1.nodes ) {
                if( pair.second.EdgeNum() > 2 )
                    return false;
            }
            return tmp_g1.NodesSize() > 2 ;
        }

        bool Simplify(GraphG1 & graph_g1 ) {
            if( ! CanSimplify(graph_g1) ) return IsLinear(graph_g1) ;
            std::vector< std::tuple< int , unsigned int > > index ;
            for( const auto & pair : graph_g1.nodes ) 
                if( pair.second.EdgeNum() > 2 ) 
                    index.push_back(std::make_tuple((int) pair.second.EdgeNum(),pair.first));
            std::sort(index.rbegin() , index.rend());
            for( const auto & tp : index )
            {
                unsigned int node_id = std::get<1>(tp);
                const auto & node = graph_g1.GetNode(node_id) ;
                if( node.EdgeNum() <= 2 ) 
                    continue ;
                std::vector<std::tuple< int , int > > edge_index ;
                for( int edge_id :  node.edge_ids ) 
                    edge_index.push_back( std::make_tuple( graph_g1.GetEdge(edge_id).weight , edge_id )) ; 
                std::sort(edge_index.rbegin(),edge_index.rend());
                if( std::get<0>(edge_index[1]) > std::get<0>(edge_index[2]) ) {
                    for( int i = 2 ; i < (int)edge_index.size() ; i ++ )
                        graph_g1.RemoveEdge(std::get<1>(edge_index[i]));
                    if( IsLinear(graph_g1) ) 
                        return true ;
                } else {
                    return false ;
                }
            }
            return false ;
        }

        void UpdateSucc( const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info ,
                const GraphG1 & graph_g1 ) 
        {

        }
        void UpdateFailed( const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info ,
                const GraphG1 & graph_g1 ) 
        {

        }
        void CorrectGraph()
        {
            mst_v1 = base_contig_sim_graph.MinTree();
            auto mst_mid = mst_v1;
            BGIQD::stLFR::ContigSimGraph::JunctionInfo junction_info 
                = mst_mid.NextJunction();
            while( junction_info.valid  )
            {
                auto graph_g1 = GetG1(junction_info) ;
                if( Simplify( graph_g1 ) ) 
                    UpdateSucc(junction_info , graph_g1);
                else 
                    UpdateFailed(junction_info , graph_g1);
                junction_info = mst_mid.NextJunction() ;
            }
            mst_v2 = base_contig_sim_graph.MinTree();
        }

        std::string log_str() const {
            return "";
        }

        std::vector<LinearOrder> GenerateLinear()
        {
            std::vector<LinearOrder> ret ;
            auto mst_linear = mst_v2 ;
            if( mst_v2.nodes.size() < 1 )
                return ret ;
            auto tip_result = BGIQD::stLFR::ContigSimGraph
                ::RemoveTip_n2(mst_linear) ;
            auto junction_result = BGIQD::stLFR::ContigSimGraph
                ::DetectJunctions(mst_linear) ;

            if( junction_result.size() > 0 )
            {
                for( auto x : junction_result )
                {
                    mst_linear.RemoveNode(x);
                }
            }
            auto splits = BGIQD::stLFR::ContigSimGraph::UnicomGraph(mst_linear);
            for( auto & pair : splits)
            {
                auto a_line = BGIQD::stLFR::ContigSimGraph::TrunkLinear(pair.second);
                if( a_line.size() > 1)
                    ret.push_back(a_line);
            }
            return ret ;
        }
    };


    BGIQD::SOAP2::FileNames fNames ;

    BGIQD::stLFR::ContigSimGraph graph;

    BGIQD::LOG::logger lger;

    std::map< BGIQD::stLFR::ContigSimGraph::NodeId , MST_correct_new>  split_graphs;

    float smallest ;

    void Init(const std::string & prefix, int min_c , float m53 , float min_js )
    {
        fNames.Init(prefix);
        smallest = min_js  ;
        min_53 = m53 ;
        min_common_barcode_type = min_c ;
        BGIQD::LOG::logfilter::singleton().get("MST", BGIQD::LOG::loglevel::INFO,lger);
    }

    void LoadContigSimGraph()
    {
        graph.use_salas = false;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.cluster("mst"));
        if( in == NULL )
            FATAL("failed to open xxx.mst.cluster for read!!! ");
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
            split_graphs[pair.first].Init(pair.second,min_common_barcode_type ,min_53);
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
        BGIQD::FREQ::Freq<int>  trunk_freq;
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

    void LoadBarcodeOnContig()
    {

    }
    //      contig_id       barcodes set
    std::map<unsigned int , std::set<int> > contig_barcodes_map;
    float min_53 ;
    int min_common_barcode_type ;
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix .\n\
                        Input xxx.cluster ; xxx.barcodeOnContig ;\n\
                        Output xxx.mintree_trunk_linear && xxx.mst_base && xxx.mst_modified ");
        DEFINE_ARG_OPTIONAL(int  , min_common_barcode_type , " min common barcode type to comfirm a road ","3");
        DEFINE_ARG_OPTIONAL(float, min_53, "min 53 ","1.5");
        DEFINE_ARG_OPTIONAL(float, min_js, "min js threshold","0.1");
    END_PARSE_ARGS
    config.Init( prefix.to_string() , min_common_barcode_type.to_int() , min_53.to_float(), min_js.to_float() );
    config.LoadContigSimGraph();
    config.LoadBarcodeOnContig();
    config.SplitGraph();
    config.CorrectGraph();
    config.GenerateLinears();
    return 0 ;
}
