#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "stLFR/CBB.h"
#include "stLFR/contigSimGraph.h"
#include "stLFR/CBB.h"

#include "soap2/fileName.h"

#include <set>
#include <vector>


struct AppConf
{
    BGIQD::SOAP2::FileNames fNames ;

    BGIQD::stLFR::ContigSimGraph graph;

    std::map<BGIQD::stLFR::ContigSimGraph::NodeId , BGIQD::stLFR::ContigSimGraph> mintrees;

    std::map<BGIQD::stLFR::ContigSimGraph::NodeId , BGIQD::stLFR::ContigSimGraph> split_graphs;

    std::map<BGIQD::stLFR::ContigSimGraph::NodeId , BGIQD::stLFR::ContigSimGraph> mintreetrunks;

    void Init(const std::string & prefix)
    {
        fNames.Init(prefix);
    }

    void LoadContigSimGraph()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.cluster());
        if( in == NULL )
            FATAL("failed to open xxx.cluster for read!!! ");
        auto parseline = [&](const std::string & line )
        {
            BGIQD::stLFR::ContigRelation tmp;
            tmp.InitFromString(line);
            for( auto x : tmp.sims )
            {
                graph.AddEdgeSim( tmp.contigId , x.first , x.second.simularity );
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
    }


    void SplitGraph()
    {
        split_graphs = graph.UnicomGraph(graph);
    }

    void GenerateMinTrees()
    {
        auto out1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mintree());
        if( out1 == NULL )
            FATAL(" failed to open xxx.mintree for write !!! ");
        for(const auto & pair : split_graphs)
        {
            auto mintree = pair.second.MinTree() ;
            mintrees[pair.first] = mintree;
            mintree.PrintAsDOT(*out1) ;
        }
        delete out1;
    }

    void GenerateMinTreeTrunks()
    {
        auto out2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mintreetrunk());
        if( out2 == NULL )
            FATAL(" failed to open xxx.mintree_trunk for write !!! ");
        for(const auto & pair : mintrees)
        {
            auto trunk = graph.TrunkFromMinTree(pair.second);
            mintreetrunks[pair.first] = trunk ;
            trunk.PrintAsDOT(*out2) ;
        }
        delete out2;
    }

    void GenerateMinTreeTrunkLinears()
    {
        int trunk_count = 1;;
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.mintreetrunklinear());
        if( out3 == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for write !!! ");
        for(const auto & pair : mintreetrunks )
        {
            (*out3)<<"---\t"<<trunk_count<<"\t---"<<std::endl;
            trunk_count++;
            auto linear = graph.TrunkLinear(pair.second);
            for(const auto x : linear )
            {
                (*out3)<<x<<std::endl;
            }
        }
        delete out3;
    }
    void test()
    {
        //BGIQD::FILES
    }
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.connInfo . Output xxx.minTree ");
    END_PARSE_ARGS

    config.Init(prefix.to_string());
    config.LoadContigSimGraph();
    config.SplitGraph();
    config.GenerateMinTrees();
    config.GenerateMinTreeTrunks();
    config.GenerateMinTreeTrunkLinears();
    return 0 ;
}
