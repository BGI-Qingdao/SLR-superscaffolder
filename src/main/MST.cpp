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

    struct MST_correct
    {

        struct simplify_log
        {

            int base_basic_node_num ;

            int base_node_num ;

            int mst_node_num ;

            BGIQD::stLFR::ContigSimGraph::TipRemoveResult tip_result;

            std::vector<BGIQD::stLFR::ContigSimGraph::NodeId> junction_result;

            void Init(int basic , int curr )
            {
                base_basic_node_num = basic  ;
                base_node_num = curr ;
                mst_node_num = 0 ;
                tip_result.Init();
                junction_result.clear();
            }
        } ;

        void Init( const  BGIQD::stLFR::ContigSimGraph & basic , int mr , float df )
        {
            base_contig_sim_graph = basic ;
            max_correct_loop_num = mr ;
            max_del_fac = df ;
            basic_num = base_contig_sim_graph.nodes.size();
        }

        BGIQD::stLFR::ContigSimGraph base_contig_sim_graph ;

        BGIQD::stLFR::ContigSimGraph maximum_span_graph ;

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

        void Simplify()
        {
            maximum_span_graph = base_contig_sim_graph.MinTree();
            curr_log.mst_node_num = maximum_span_graph.nodes.size() ;
            if ( curr_log.mst_node_num > 3 )
            {
                curr_log.tip_result = BGIQD::stLFR::ContigSimGraph
                    ::RemoveTip_n2(maximum_span_graph) ;
                curr_log.junction_result = BGIQD::stLFR::ContigSimGraph
                    ::DetectJunctions(maximum_span_graph) ;

                if( curr_log.junction_result.size() > 0 )
                {
                    for( auto x : curr_log.junction_result )
                        base_contig_sim_graph.RemoveNode(x);
                }
            }
        }

        bool IsClean() const 
        {
            if ( logs.empty() )
                return false ;

            const auto & item = logs.top() ;

            return item.base_node_num < 4
                || ( int(logs.size()) >= max_correct_loop_num )
                || ( item.base_node_num * 10 / item.base_basic_node_num <= 9 )
                || ( item.tip_result.tip_num == 0 && item.junction_result.size() == 0 ) ;
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
                    maximum_span_graph.RemoveNode(x);
            }
            auto splits = BGIQD::stLFR::ContigSimGraph::UnicomGraph(maximum_span_graph);
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

    std::map< BGIQD::stLFR::ContigSimGraph::NodeId , MST_correct>  split_graphs;

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
        graph.use_salas = false;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.cluster());
        std::vector<std::tuple<float,unsigned int , unsigned int> >simcache;
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
                    simcache.push_back(std::make_tuple( x.second.simularity,tmp.contigId , x.first ));
                }
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
        lger<<BGIQD::LOG::lstart() << "load contig sim graph done "<<BGIQD::LOG::lend() ;
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
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.minTree ");
        DEFINE_ARG_OPTIONAL(float , threshold, " threshold of simularity","0.1");
        DEFINE_ARG_OPTIONAL(float , del_fac, " threshold of del junction node","0.9");
        DEFINE_ARG_OPTIONAL(int, del_round, "maximum del round ","3");
    END_PARSE_ARGS

    config.Init(prefix.to_string(),threshold.to_float());
    config.max_del_fac = del_fac.to_float();
    config.max_correct_loop_num = del_round.to_int() ;
    config.LoadContigSimGraph();
    config.SplitGraph();
    config.CorrectGraph();
    config.GenerateLinears();
    return 0 ;
}
