#ifndef __STLFR_CONTIGCLUSTER_H__
#define __STLFR_CONTIGCLUSTER_H__

#include <map>

namespace BGIQD{
    namespace stLFR{


        struct ContigCluster
        {
            typedef std::map< int , std::map< int , float > >  cluters;
            cluters connections;
            int clusterNum;
            void loadCluster(const std::string & file);
        };


    }
}
#endif //__STLFR_CONTIGCLUSTER_H__
