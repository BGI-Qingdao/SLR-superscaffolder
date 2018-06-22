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
    mintree.PrintAsDOT() ;
    return 0 ;
}
