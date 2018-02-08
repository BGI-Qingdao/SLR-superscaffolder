#ifndef __SOAP2_LOADGRAPH_H__
#define __SOAP2_LOADGRAPH_H__

#include <string>
#include "soap2/contigGraph.h"

namespace BGIQD {
namespace SOAP2 {

struct GlobalConfig
{
    Edge * edge_array;
    Arc * arc_array;
    KeyEdge * key_array;
    std::string updateEdge;
    std::string arc;
    std::string cluster;
    unsigned int contigTotalNum;
    unsigned int contigNum;
    unsigned int clusterNum;
    long long arcNum;
    long long connectionNum;

    std::map<unsigned int , unsigned int> connections;
};

//void loadContigIndex(Edge * array , GlobalConfig & config);
//void loadContig(Edge * array);
//
// 
void loadUpdateEdge(Edge * array , GlobalConfig & config);
void loadArc( Edge * array , GlobalConfig & config);
void loadCluster(KeyEdge * array , GlobalConfig & config);
void buildConnection(KeyEdge * array , Edge * e_array);

}
}
#endif //__SOAP2_LOADGRAPH_H__
