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
        bool ContainCircle( const GraphG1 & graph)
        {
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
        BGIQD::stLFR::ContigSimGraph maximum_span_graph ;
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

        void CorrectGraph()
        {
            maximum_span_graph = base_contig_sim_graph.MinTree();
            BGIQD::stLFR::ContigSimGraph::JunctionInfo junction_info 
                = base_contig_sim_graph.NextJunction();
            while( junction_info.valid  )
            {
                auto graph_g1 = GetG1(junction_info) ;

                junction_info = base_contig_sim_graph.NextJunction() ;
            }
            maximum_span_graph = base_contig_sim_graph.MinTree();
        }

        std::string log_str() const {
            return "";
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
    //config.GenerateLinears();
    return 0 ;
}
