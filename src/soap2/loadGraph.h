#ifndef __SOAP2_LOADGRAPH_H__
#define __SOAP2_LOADGRAPH_H__

#include <string>
#include "soap2/contigGraph.h"
#include <set>
#include <mutex>

namespace BGIQD {
    namespace SOAP2 {

        struct ContigRoad
        {
            bool headin;
            bool tailin;
            int length;
            bool downstream;
            std::vector<unsigned int> contig;
            std::vector<unsigned int> real_contig;
        };


        struct GlobalConfig
        {
            // graph info
            GraphEA  graph_ea;
            Edge * edge_array;
            //Arc * arc_array;
            //long long arcNum;
            //unsigned int contigTotalNum;

            // base config
            int K;
            std::string updateEdge;
            std::string arc;
            std::string cluster;
            std::string contigroad;
            // cluster(key) data
            unsigned int clusterNum;
            std::set<unsigned int> keys;
            std::map<unsigned int , std::map<unsigned int,float > > connections;
            long long connectionNum;

            // edge id --> key id
            std::map<unsigned int , unsigned int> key_map;
            KeyEdge * key_array;
            std::mutex * key_mutex;

            // super conitgs
            std::mutex contig_mutex;
            std::vector<ContigRoad> contigs ;//std::vector<unsigned int> > contigs;
        };

        void loadUpdateEdge( GlobalConfig & config);
        void loadArc( GlobalConfig & config);
        void loadCluster(GlobalConfig & config);
        //void buildConnection(GlobalConfig&);
        //void LinearConnection(GlobalConfig &);
    }
}
#endif //__SOAP2_LOADGRAPH_H__
