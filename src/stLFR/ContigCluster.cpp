#include "stLFR/ContigCluster.h"
#include "common/files/file_reader.h"
#include <sstream>

namespace BGIQD{
    namespace stLFR{

        void ContigCluster::loadCluster( const std::string & file)
        {
            std::string line;
            unsigned int contigId;
            unsigned int to;
            float cov;
            clusterNum = 0;
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            // load connection 
            while(!std::getline(*in,line).eof())
            {
                std::istringstream ist(line);
                ist>>contigId;

                while(! ist.eof() )
                {
                    ist>>to>>cov;
                    connections[contigId][to] = cov;
                    connections[to][contigId] = cov ;
                }
            }
            delete in ;
            clusterNum = connections.size();
        }
    }
}
