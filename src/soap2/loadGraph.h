#ifndef __SOAP2_LOADGRAPH_H__
#define __SOAP2_LOADGRAPH_H__

#include <string>
#include "soap2/contigGraph.h"
#include <set>

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
    unsigned int clusterNum;
    long long arcNum;
    long long connectionNum;

    std::map<unsigned int , std::map<unsigned int,float > > connections;
    std::set<unsigned int> keys;
};

//void loadContigIndex(Edge * array , GlobalConfig & config);
//void loadContig(Edge * array);
//
// 
void loadUpdateEdge( GlobalConfig & config);
void loadArc( GlobalConfig & config);
void loadCluster(GlobalConfig & config);
void buildConnection(GlobalConfig&);

}
}
#endif //__SOAP2_LOADGRAPH_H__
