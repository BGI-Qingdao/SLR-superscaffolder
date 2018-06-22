#include "common/args/argsparser.h"
#include "stLFR/CBB.h"
#include "soap2/fileName.h"
#include "stLFR/contigSimGraph.h"

#include <set>
#include <vector>


struct AppConf
{
}config;


int main(int argc , char **argv )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.connInfo . Output xxx.minTree ");
    END_PARSE_ARGS

    BGIQD::stLFR::ContigSimGraph graph;
    auto mintree = graph.MinTree() ;
    mintree.PrintAsDOT() ;
    return 0 ;
}
