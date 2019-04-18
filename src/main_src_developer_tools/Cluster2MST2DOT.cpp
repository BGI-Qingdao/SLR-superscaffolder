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

    BGIQD::SOAP2::FileNames fNames ;
    BGIQD::stLFR::ContigSimGraph graph;
    BGIQD::LOG::logger lger;

    float smallest ;

    void Init(const std::string & prefix, float f)
    {
        fNames.Init(prefix);
        smallest = f ;
        BGIQD::LOG::logfilter::singleton().get("MST2DOT", BGIQD::LOG::loglevel::INFO,lger);
    }

    void LoadContigSimGraph()
    {
        graph.use_salas = false;
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

    void GenMST()
    {
        auto mst = graph.MinTree();
        mst.PrintAsDOT(std::cout);
    }

}config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.minTree ");
        DEFINE_ARG_OPTIONAL(float , threshold, " threshold of simularity","0.1");
    END_PARSE_ARGS

    config.Init(prefix.to_string(),threshold.to_float());
    config.LoadContigSimGraph();
    config.GenMST();
    return 0 ;
}
