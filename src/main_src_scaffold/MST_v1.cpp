#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"
#include "common/stl/setHelper.h"
#include "stLFR/CBB.h"
#include "stLFR/contigSimGraph.h"
#include "stLFR/CBB.h"

#include "soap2/fileName.h"
#include "soap2/contigIndex.h"
#include <set>
#include <vector>
#include <stack>
#include <tuple>
#include <sstream>

BGIQD::FREQ::Freq<std::string>  Simplify_freq;
BGIQD::FREQ::Freq<std::string>  Failed_reason_freq;
BGIQD::FREQ::Freq<std::string>  Cant_reason_freq;
BGIQD::FREQ::Freq<std::string>  Road_result_freq;
std::set<unsigned int> masked_nodes ;
std::string log_str() {
	return "\nSimplify result:\n"+ Simplify_freq.ToString() + "\n" + 
		"Simplify failed reason:\n"+ Failed_reason_freq.ToString() + "\n" + 
		"Can't Simplify reason:\n"+ Cant_reason_freq.ToString() + "\n" + 
        "Check path result : \n" +Road_result_freq.ToString() + "\n";
}
std::map<unsigned int , float> lf_log ;
std::map<unsigned int , float> rf_log ;
std::map<unsigned int , float> mf_log ;
void PrintSingloten( const std::string & file) {
    auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
    if(out==NULL)
        FATAL(" can't open singleton file for write !");
    for( const auto & pair : lf_log) 
        (*out)<<"L\t"<<pair.first<<'\t'<<pair.second<<'\n';
    for( const auto & pair : mf_log) 
        (*out)<<"M\t"<<pair.first<<'\t'<<pair.second<<'\n';
    for( const auto & pair : rf_log) 
        (*out)<<"R\t"<<pair.first<<'\t'<<pair.second<<'\n';
    delete out ;
}
std::map<unsigned int ,BGIQD::SOAP2::ContigIndex> seeds;
struct AppConf
{
    typedef std::vector<BGIQD::stLFR::ContigSimGraph::NodeId> LinearOrder;

    struct MST_correct_new 
    {
        float del_fac ;
        int del_round ;

        const std::map<unsigned int , std::set<int> > * contig_barcodes_map_ptr ;
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
        struct BarcodeTypeInfo : public OldRoad {
            int left_only ;
            int mid_only ;
            int right_only ;
            int common ;
        };

        // Done
        BarcodeTypeInfo GetBarcodeDetail(const OldRoad & old_road ) const {
            BarcodeTypeInfo ret ;
            const auto & left_set = contig_barcodes_map_ptr->at(old_road.left) ; 
            const auto & mid_set = contig_barcodes_map_ptr->at(old_road.mid) ; 
            const auto & right_set = contig_barcodes_map_ptr->at(old_road.right) ; 
	
            ret.left = old_road.left ;
            ret.mid = old_road.mid;
            ret.right = old_road.right;
            auto comm1 = BGIQD::STL::set_common(left_set,mid_set);
            auto comm = BGIQD::STL::set_common(comm1,right_set);
            ret.common = comm.size() ;

            auto left_diff1 = BGIQD::STL::set_diff_in_s1(left_set,mid_set);
            auto left_diff = BGIQD::STL::set_diff_in_s1(left_diff1,right_set);
            ret.left_only = left_diff.size();

            auto mid_diff1 = BGIQD::STL::set_diff_in_s1(mid_set,left_set);
            auto mid_diff = BGIQD::STL::set_diff_in_s1(mid_diff1,right_set);
            ret.mid_only = mid_diff.size();

            auto right_diff1 = BGIQD::STL::set_diff_in_s1(right_set,left_set);
            auto right_diff = BGIQD::STL::set_diff_in_s1(right_diff1,mid_set);
            ret.right_only = right_diff.size();

            return ret ;
        }

        // Done
        static std::vector<OldRoad> GetCandidateRoads(  const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info ) {
            std::vector<OldRoad> ret ;
            for( int i = 0 ; i < (int)junction_info.neibs.size() ; i ++ )
            {
                for( int j = i+1 ; j < (int)junction_info.neibs.size() ; j ++ )
                    ret.emplace_back( OldRoad{ junction_info.neibs.at(i) 
                            , junction_info.junction_id
                            , junction_info.neibs.at(j) } );

            }
            return ret ;
        }
        //Done
        NewRoad CheckRoad( const OldRoad & old_road ) {
            NewRoad ret ;
            ret.valid = false ;
            auto barcode_type_info = GetBarcodeDetail(old_road);
            if( barcode_type_info.common < min_common_barcode_type ) {
                Road_result_freq.Touch("No enough common");
                return ret ;
            }
            std::vector< std::tuple< float ,unsigned  int > > index ;
            float lf , mf , rf ;
            lf = float(barcode_type_info.left_only) / float(seeds.at(barcode_type_info.left).length); 
            rf = float(barcode_type_info.right_only) / float(seeds.at(barcode_type_info.right).length); 
            mf = float(barcode_type_info.mid_only) / float(seeds.at(barcode_type_info.mid).length); 
            index.push_back(std::make_tuple( lf , barcode_type_info.left) );
            index.push_back(std::make_tuple( mf , barcode_type_info.mid ) );
            index.push_back(std::make_tuple( rf , barcode_type_info.right) );
            std::sort(index.rbegin() , index.rend());
            float tmp_bn1 , tmp_bn2 ,tmp_bn3 ;
            unsigned int tmp_ci1 , tmp_ci2 , tmp_ci3 ;
            std::tie(tmp_bn3 , tmp_ci3) = index[2] ;
            std::tie(tmp_bn2 , tmp_ci2) = index[1] ;
            std::tie(tmp_bn1 , tmp_ci1) = index[0] ;
            if( tmp_bn3 <= 0 ) {
                Road_result_freq.Touch("No monoply barcode");
                return ret ;
            }
            if( float(tmp_bn2) / float(tmp_bn3) < min_53 ) {
                Road_result_freq.Touch("No enough differece of 2&3");
                return ret ;
            }
            Road_result_freq.Touch("Succ");
            ret.valid = true ;
            ret.left = tmp_ci1;
            ret.mid  = tmp_ci3; // put the smallest in mid
            ret.right = tmp_ci2 ;
            return ret ;
        }
        //Done
        static bool ContainCircle( const GraphG1 & graph){
            auto splits = graph.UnicomGraph() ;
            for( const auto & pair : splits ){
                const auto & a_graph = pair.second ;
                if( a_graph.NodesSize() != a_graph.EdgesSize() +1 )
                    return true;
            }
            return false ;
        }

        //Done
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
            //if( ! ContainCircle(ret ) ) {
            //for( auto neib : junction_info.neibs )
            //    ret.Add1Edge(neib,junction_info.junction_id);
            //}
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
        int round ;
        int del_node ;
        struct MSTLinear{
            unsigned int mid ;
            char is_smallest ;
            int common_barcode_num ;
            float f23 ;
        };
        std::vector<MSTLinear> GetMSTLinearInfo( const std::map<unsigned int ,BGIQD::SOAP2::ContigIndex> & cls) {
            mst_v1 = base_contig_sim_graph.MinTree();
            std::vector<MSTLinear> ret;
            for( const auto & pair : mst_v1.nodes ) {
                const auto & node = pair.second ;
                if ( node.EdgeNum() != 2 ) continue ;
                std::vector<unsigned int > lr ;
                for( int edge_id : node.edge_ids ) {
                    const auto & edge = mst_v1.GetEdge(edge_id);
                    const auto & anode = mst_v1.GetNode( edge.OppoNode(node.id) );
                    if( anode.EdgeNum() == 2 ) lr.push_back( anode.id ) ;
                }
                if( lr.size() != 2 ) continue ;
                OldRoad tmp ;
                tmp.left = lr[0] ; tmp.mid = node.id ; tmp.right = lr[1] ;
                auto detail = GetBarcodeDetail(tmp);
                float f23 = 0.0f ;
                float lr_min = 0;
                float lf , rf , mf ;
                lf = float(detail.left_only) / float(cls.at(detail.left).length); 
                rf = float(detail.right_only) / float(cls.at(detail.right).length); 
                mf = float(detail.mid_only) / float(cls.at(detail.mid).length); 
                lf_log[detail.left] = lf ;
                rf_log[detail.right] = rf ;
                mf_log[detail.mid] = mf ;
                if( lf < rf ) 
                    lr_min = lf ;
                else
                    lr_min = rf;
                if (  mf != 0 )
                    f23 = lr_min / mf ;
                ret.push_back( { node.id ,(mf<=lr_min ? 'Y':'N'), detail.common , f23 } );
            }
            return ret ;
        }
        // Done
        void Init(const BGIQD::stLFR::ContigSimGraph & base , int mc , float m53 ,
                const std::map<unsigned int , std::set<int> > & map
                ) {
            base_contig_sim_graph = base ;
            min_53 = m53 ;
            min_common_barcode_type = mc ;
            base_num = base.NodesSize() ; 
            contig_barcodes_map_ptr = &map;
        }
        //Done
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
            } else if ( splits.size() ==1 ) {
                group_num = 1 ;
                the_one = &(splits.begin()->second);
            } else {
                assert(0);
            }
            assert(group_num > 0 && the_one != NULL );
            return std::make_pair( group_num , *the_one ) ;
        }
        //Done
        static bool CanSimplify( const  GraphG1 & graph_g1 ) {
            if( graph_g1.EdgesSize() < 2 ) {
                  Cant_reason_freq.Touch("Edge<2");
                  return false ;
            }
            auto gret = GetA_group(graph_g1) ;
            if( gret.first > 1  ) {
               Cant_reason_freq.Touch("SubGraph>1");
               return false ;
            }
            const auto & tmp_g1 = gret.second ;
            for( const auto & pair : tmp_g1.nodes ) {
                if( pair.second.EdgeNum() > 2 )
                    return true ;
            }
            Cant_reason_freq.Touch("NoJunction");
            return false;
        }
        //Done
        static bool IsLinear( const  GraphG1 & graph_g1) {
            if( graph_g1.EdgesSize() < 2 ) return false ;
            auto gret = GetA_group(graph_g1) ;
            if( gret.first > 1  ) return false ;
            const auto & tmp_g1 = gret.second ;
            for( const auto & pair : tmp_g1.nodes ) {
                if( pair.second.EdgeNum() > 2 )
                    return false;
            }
            return tmp_g1.NodesSize() > 2 ;
        }
        //Done
        static bool Simplify(GraphG1 & graph_g1 ) {
            if( ! CanSimplify(graph_g1) ) {
                Failed_reason_freq.Touch("Can't Simplify");
                return IsLinear(graph_g1) ;
            }
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
                for( int edge_id :  node.edge_ids ) {
                    auto & edge = graph_g1.GetEdge(edge_id);
                    edge_index.push_back( std::make_tuple( edge.weight , edge_id )) ; 
                }
                std::sort(edge_index.rbegin(),edge_index.rend());
                if( std::get<0>(edge_index[1]) > std::get<0>(edge_index[2]) ) {
                    for( int i = 2 ; i < (int)edge_index.size() ; i ++ )
                        graph_g1.RemoveEdge(std::get<1>(edge_index[i]));
                    if( IsLinear(graph_g1) ) 
                        return true ;
                } else {
                    Failed_reason_freq.Touch("Edge same weight");
                    return false ;
                }
            }
            Failed_reason_freq.Touch("Final Failed");
            return false ;
        }
        //Done
        static std::vector<std::pair< int , int > > GetLinearFromG1(const GraphG1 & g1 ){
            std::vector<std::pair<int ,int> > ret ;
            auto gret = GetA_group(g1);
            assert( gret.first == 1  );
            const auto & a_graph = gret.second ;
            for( const auto & edge :  a_graph.edges) 
                if( edge.IsValid())
                    ret.push_back( std::make_pair( edge.from, edge.to ) );
            return ret ;
        }
        //Done
        static void UpdateSucc(  const GraphG1 & graph_g1 ,
                BGIQD::stLFR::ContigSimGraph & mst ,
                BGIQD::stLFR::ContigSimGraph & contig_sim
                )
        {
            auto new_edges = GetLinearFromG1(graph_g1) ;
            for( const auto & pair :new_edges ) 
                contig_sim.AddEdgeSim(pair.first,pair.second,1);
            for( const auto & pair :new_edges ) 
                mst.AddEdgeSim(pair.first,pair.second,1);
        }
        //Done
        void DeleteJunctions( const BGIQD::stLFR::ContigSimGraph::JunctionInfo & junction_info ,
                BGIQD::stLFR::ContigSimGraph & mst ,
                BGIQD::stLFR::ContigSimGraph & contig_sim
                ) 
        {
            contig_sim.RemoveNode(junction_info.junction_id);
            mst.RemoveNode(junction_info.junction_id);
        }
        static void MaxMSTEdge(
                BGIQD::stLFR::ContigSimGraph & mst ,
                BGIQD::stLFR::ContigSimGraph & contig_sim
                )
        {
            for( const auto & edge : mst.edges ) 
                if( edge.IsValid())
                    contig_sim.UpdateEdge(edge.from ,edge.to , 1.0 );
        }
        bool IsComplexJunction( 
                const BGIQD::stLFR::ContigSimGraph & mst  ,
                const BGIQD::stLFR::ContigSimGraph::JunctionInfo & i ) const {
            const auto & node = mst.GetNode(i.junction_id) ;
            return node.EdgeNum() > 4 ;
        }
        //Done
        void CorrectGraph()
        {
            del_node = 0 ;
            round = 0 ;
            mst_v1 = base_contig_sim_graph.MinTree();
            while( ( 1.0f - float(del_node)/ float(base_num) ) > del_fac && round<= del_round  ) {
                round ++ ;
                auto mst_mid = base_contig_sim_graph.MinTree();
                BGIQD::stLFR::ContigSimGraph::RemoveTip_n2(mst_mid) ;
                BGIQD::stLFR::ContigSimGraph::JunctionInfo junction_info 
                    = mst_mid.NextJunction();
                bool change = false ;
                while( junction_info.valid  )
                {
                    change = true ;
                    if( IsComplexJunction( mst_mid ,junction_info ) ) {
                        DeleteJunctions( junction_info , mst_mid , base_contig_sim_graph );
                        del_node ++ ;
                        masked_nodes.insert(junction_info.junction_id);
                        junction_info = mst_mid.NextJunction();
                        continue ;
                    }
                    DeleteJunctions( junction_info , mst_mid , base_contig_sim_graph );
                    auto graph_g1 = GetG1(junction_info) ;
                    if( Simplify( graph_g1 ) ) {
                        UpdateSucc( graph_g1 , mst_mid , base_contig_sim_graph );
                        Simplify_freq.Touch("Simplify succ");
                    } else {
                        del_node ++ ;
                        Simplify_freq.Touch("Simplify failed");
                        masked_nodes.insert(junction_info.junction_id);
                    }
                    junction_info = mst_mid.NextJunction();
                }
                if( ! change ) break ;
            }
            mst_v2 = base_contig_sim_graph.MinTree();
        }
        struct TipRmLinear {
            bool is_tip ;
            std::vector<unsigned int> linear_nib ;
            std::vector<unsigned int> tip_nib ;
            std::vector<unsigned int> junc_nib ;
        };
        static TipRmLinear TipCheck(
                const BGIQD::stLFR::ContigSimGraph & mst  ,
                const BGIQD::stLFR::ContigSimGraph::JunctionInfo & i ) {
            TipRmLinear ret ;
            const auto & node = mst.GetNode(i.junction_id);
            for( int edge_id : node.edge_ids ) {
                const auto & nib = mst.GetNode( mst.GetEdge(edge_id).OppoNode(node.id) );
                if( nib.EdgeNum() == 1 ) ret.tip_nib.push_back(nib.id);
                if( nib.EdgeNum() == 2 ) ret.linear_nib.push_back(nib.id);
                if( nib.EdgeNum() >2 )   ret.junc_nib.push_back(nib.id);
            }
            if( ret.junc_nib.size() == 0 && ret.linear_nib.size() == 2 ) ret.is_tip = true ;
            else ret.is_tip = false ;
            return ret ;
        }

        //Done
        static void UpdateTip( 
                unsigned int junction_id ,
                const TipRmLinear & tip,
                BGIQD::stLFR::ContigSimGraph & mst
                )
        {
            assert( tip.linear_nib.size() == 2 ) ;
            for( const auto & nib : tip.linear_nib ) 
                mst.AddEdgeSim(nib,junction_id,1.0f);
        }


        //Done
        std::vector<LinearOrder> GenerateLinear()
        {
            std::vector<LinearOrder> ret ;
            auto mst_linear = mst_v2 ;
            if( mst_v2.nodes.size() < 1 )
                return ret ;

            BGIQD::stLFR::ContigSimGraph::JunctionInfo junction_info 
                = mst_linear.NextJunction();
            while( junction_info.valid  )
            {
                auto ret = TipCheck( mst_linear , junction_info );
                if( ! ret.is_tip && IsComplexJunction( mst_linear,junction_info ) ) {
                    DeleteJunctions( junction_info , mst_linear, base_contig_sim_graph );
                    masked_nodes.insert(junction_info.junction_id);
                    junction_info = mst_linear.NextJunction();
                    continue ;
                }
                DeleteJunctions( junction_info , mst_linear, base_contig_sim_graph );
                auto graph_g1 = GetG1(junction_info) ;
                if( Simplify( graph_g1 ) ) {
                    UpdateSucc( graph_g1 , mst_linear , base_contig_sim_graph );
                    Simplify_freq.Touch("Simplify succ");
                } else {
                    if( ret.is_tip ) {
                        UpdateTip(junction_info.junction_id , ret , mst_linear);
                    }else {
                        Simplify_freq.Touch("Simplify failed");
                        masked_nodes.insert(junction_info.junction_id);
                    }
                }
                junction_info = mst_linear.NextJunction();
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

    //Done
    void Init(const std::string & prefix, int min_c , float m53 , float min_js )
    {
        fNames.Init(prefix);
        smallest = min_js  ;
        min_53 = m53 ;
        min_common_barcode_type = min_c ;
        BGIQD::LOG::logfilter::singleton().get("MST_v1", BGIQD::LOG::loglevel::INFO,lger);
    }

    // Done
    void LoadContigSimGraph()
    {
        graph.use_salas = false;
        lger<<BGIQD::LOG::lstart() << "load contig sim graph start ..."<<BGIQD::LOG::lend() ;
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
        lger<<BGIQD::LOG::lstart() << "contig-sim graph nodes : "<<graph.NodesSize()<<BGIQD::LOG::lend() ;
        lger<<BGIQD::LOG::lstart() << "contig-sim graph edges : "<<graph.EdgesSize()<<BGIQD::LOG::lend() ;
    }

    // Done
    void SplitGraph()
    {
        auto splits = graph.UnicomGraph(graph);
        for( const auto & pair : splits )
        {
            split_graphs[pair.first].Init(pair.second,min_common_barcode_type ,min_53,contig_barcodes_map);
            split_graphs[pair.first].del_fac = del_fac ;
            split_graphs[pair.first].del_round = del_round ;
        }
        lger<<BGIQD::LOG::lstart() << "split contig sim graph into "<<split_graphs.size()<<" sub graph"<<BGIQD::LOG::lend() ;
    }

    //Done
    void CorrectGraph()
    {
        for( auto & pair : split_graphs)
        {
            pair.second.CorrectGraph();
        }
        lger<<BGIQD::LOG::lstart() <<"\n" <<log_str()<<BGIQD::LOG::lend() ;
    }
    //Done
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

    // Done
    void LoadBarcodeOnContig()
    {
        lger<<BGIQD::LOG::lstart() << "load barcodeOnContig start ..."<<BGIQD::LOG::lend() ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.BarcodeOnContig());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnContig to read !");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            if( ! graph.HasNode( tmp.contig_id ) )
                continue ;
            for( const auto & pair : tmp.barcodesOnPos)
                for( int barcode : pair.second)
                    contig_barcodes_map[tmp.contig_id].insert(barcode);
        }
        delete in;
        lger<<BGIQD::LOG::lstart() << "load barcodeOnContig done"<<BGIQD::LOG::lend() ;
    }
    //      contig_id       barcodes set
    std::map<unsigned int , std::set<int> > contig_barcodes_map;
    float min_53 ;
    int min_common_barcode_type ;

    void PrintBaiscMST(const std::string & output)
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output);
        if(out==NULL)
            FATAL(" can't open output file for write !");
        (*out) << "graph {"<<'\n';
        for(const auto & pair : split_graphs)
        {
            pair.second.mst_v1.PrintDOTEdges(*out);
        }
        (*out) << "}\n";
        delete out ;
    }

    void PrintFinalMST(const std::string & output)
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output);
        if(out==NULL)
            FATAL(" can't open output file for write !");
        (*out) << "graph {"<<'\n';
        for(const auto & pair : split_graphs)
        {
            pair.second.mst_v2.PrintDOTEdges(*out);
        }
        (*out) << "}\n";
        delete out ;
    }

    void PrintJunctionNodes()
    {
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mst_error());
        if( out3 == NULL )
            FATAL(" failed to open xxx.mst_error for write !!! ");
        for( auto & i: masked_nodes)
        {
            (*out3)<<i<<'\n';
        }
        delete out3;
    }

    void PrintLinearInfo(const std::string & file ) {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
        if(out==NULL)
            FATAL(" can't open output file for write !");
        for( auto & split : split_graphs ) {
            auto ret = split.second.GetMSTLinearInfo(seeds) ;
            for( const auto & i : ret ) 
                (*out)<<i.mid<<'\t'<<i.is_smallest<<'\t'<<i.common_barcode_num<<'\t'<<i.f23<<'\n';
        }
        delete out ;
    }
    void LoadContigIndex(){
        BGIQD::LOG::timer t(lger,"LoadContigIndex");
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.ContigIndex()) ;
        if( in == NULL )
            FATAL( "open .ContigIndex file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() ){
            BGIQD::SOAP2::ContigIndex tmp;
            tmp.InitFromString(line);
            seeds[tmp.contig] = tmp;
        }
        delete in ;
    };
    float del_fac ;
    int del_round ;
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix .\n\
                        Input xxx.cluster ; xxx.barcodeOnContig ;\n\
                        Output      xxx.mintree_trunk_linear\n\
                                    xxx.mst_base \n\
                                    xxx.mst_modified \n\
                                    xxx.mst_error");
        DEFINE_ARG_OPTIONAL(int  , min_common_barcode_type , " min common barcode type to comfirm a road ","10");
        DEFINE_ARG_OPTIONAL(float, min_53, "min 53 ","2");
        DEFINE_ARG_OPTIONAL(float, min_js, "min js threshold","0.1");
        DEFINE_ARG_OPTIONAL(bool,  debug_linear, "print linear info and exit","No");
        DEFINE_ARG_OPTIONAL(float , del_fac, " threshold of del junction node","0.9");
        DEFINE_ARG_OPTIONAL(int, del_round, "maximum del round ","1000");
    END_PARSE_ARGS
    config.Init( prefix.to_string() , min_common_barcode_type.to_int() , min_53.to_float(), min_js.to_float() );
    config.del_fac = del_fac.to_float() ;
    config.del_round = del_round.to_int() ;
    config.LoadContigIndex();
    config.LoadContigSimGraph();
    config.LoadBarcodeOnContig();
    config.SplitGraph();
    if( debug_linear.to_bool() ) {
        config.PrintLinearInfo(prefix.to_string()+".mst_linear_info" );
        PrintSingloten(prefix.to_string()+".singloten_factor");
        return 0 ;
    }
    config.CorrectGraph();
    config.GenerateLinears();
    config.PrintBaiscMST(prefix.to_string()+".mst_v1");
    config.PrintFinalMST(prefix.to_string()+".mst_v2");
    config.PrintJunctionNodes();
    return 0 ;
}
