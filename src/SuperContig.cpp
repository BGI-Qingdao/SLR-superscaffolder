#include "soap2/contigGraph.h"
#include "soap2/loadGraph.h"
//#include ""

int main()
{
    BGIQD::SOAP2::GlobalConfig config;

    BGIQD::SOAP2::loadUpdateEdge(config);
    BGIQD::SOAP2::loadArc(config);
    BGIQD::SOAP2::loadCluster(config);
    BGIQD::SOAP2::buildConnection(config);
    return 0;
}
