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

    struct MST_correct{
        void Init( const  BGIQD::stLFR::ContigSimGraph & basic  )
        {
            base_contig_sim_graph = basic ;
            basic_num = base_contig_sim_graph.nodes.size();
        }

        BGIQD::stLFR::ContigSimGraph base_contig_sim_graph ;

        BGIQD::stLFR::ContigSimGraph maximum_span_graph ;

        int basic_num ;

        void MinTree( )
        {
            maximum_span_graph = base_contig_sim_graph.MinTree();
        }

    };

    BGIQD::stLFR::ContigSimGraph graph;

    std::map< BGIQD::stLFR::ContigSimGraph::NodeId , MST_correct>  split_graphs;

    float smallest ;
    std::string cluster ;
    void Init(const std::string & c, float f)
    {
        cluster = c ;
        smallest = f ;
    }

    void LoadContigSimGraph()
    {
        graph.use_salas = false;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(cluster);
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
    }

    void SplitGraph()
    {
        auto splits = graph.UnicomGraph(graph);
        for(const auto & pair : splits)
        {
            split_graphs[pair.first].Init(pair.second);
        }
    }

    void MinTree()
    {
        for( auto & pair : split_graphs)
        {
            pair.second.MinTree();
        }
    }

    void PrintMinTree(const std::string & output)
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output);
        if(out)
            FATAL(" can't open output file for write !");
        (*out) << "graph {"<<'\n';
        for( auto & pair : split_graphs)
        {
            pair.second.maximum_span_graph.PrintDOTEdges(*out);
        }
        (*out) << "}\n";
        delete out ;
    }

}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , cluster, "Input xxx.mst.cluster .");
        DEFINE_ARG_REQUIRED(std::string , output , "Output out.mintree");
        DEFINE_ARG_OPTIONAL(float , threshold, " threshold of simularity","0.1");
    END_PARSE_ARGS

    config.Init(cluster.to_string(),threshold.to_float());
    config.LoadContigSimGraph();
    config.SplitGraph();
    config.MinTree();
    config.PrintMinTree(output.to_string());
    return 0 ;
}
