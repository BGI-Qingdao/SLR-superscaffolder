#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"

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
}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.connInfo . Output xxx.minTree ");
    END_PARSE_ARGS
    config.Init(prefix.to_string());
    config.LoadContigSimGraph();

    auto mintree = config.graph.MinTree() ;
    auto out1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(config.fNames.mintree());
    if( out1 == NULL )
        FATAL(" failed to open xxx.mintree for write !!! ");
    mintree.PrintAsDOT(*out1) ;
    delete out1;

    auto trunk = config.graph.TrunkFromMinTree(mintree);
    auto out2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(config.fNames.mintreetrunk());
    if( out2 == NULL )
        FATAL(" failed to open xxx.mintree_trunk for write !!! ");
    trunk.PrintAsDOT(*out2) ;
    delete out2;

    auto linear = config.graph.TrunkLinear(trunk);
    auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(config.fNames.mintreetrunklinear());
    if( out3 == NULL )
        FATAL(" failed to open xxx.mintree_trunk_linear for write !!! ");
    for(const auto x : linear )
    {
        (*out3)<<x<<std::endl;
    }
    delete out3;
    return 0 ;
}
