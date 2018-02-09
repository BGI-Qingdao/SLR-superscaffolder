#include "soap2/loadGraph.h"
#include <cassert>
#include <common/files/file_reader.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include "soap2/contigGraph.h"
namespace BGIQD {
namespace SOAP2 {

void loadUpdateEdge( GlobalConfig & config)
{
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(config.updateEdge);
    std::string line;
    std::getline(*in,line);
    sscanf(line.c_str(),"EDGES %u",&config.contigTotalNum);
    config.edge_array = static_cast<Edge*>( calloc(sizeof(Edge) , config.contigTotalNum + 1 ));
    int length;
    int bal;
    int cov;
    unsigned int index = 1 ;
    while(!std::getline(*in,line).eof())
    {
        sscanf(line.c_str(),">length %d,%d,%d",&length,&bal,&cov);
        config.edge_array[index].id = index ;
        config.edge_array[index].bal_id = index+bal;
        config.edge_array[index].cov = cov;
        config.edge_array[index].length = length;
        index ++ ;
    }
    assert( index == config.contigTotalNum +1 );
    return ;
}

void loadArc(GlobalConfig & config)
{
    std::string line;
    unsigned int contigId;
    unsigned int to;
    int cov;
    config.arcNum = 0;
    // Counting arcs
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(config.updateEdge);
    while(!std::getline(*in,line).eof())
    {
        std::istringstream ist(line);
        ist>>contigId;
        while(! ist.eof() )
        {
            ist>>to>>cov;
            config.arcNum ++ ;
        }
    }
    delete in ;
    config.arc_array =static_cast<Arc*>( calloc(sizeof(Arc),config.arcNum));
    long long index = 1 ;
    auto in1 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(config.updateEdge);
    while(!std::getline(*in1,line).eof())
    {
        std::istringstream ist(line);
        ist>>contigId;
        while(! ist.eof() )
        {
            ist>>to>>cov;
            config.arc_array[index].to = to;
            config.arc_array[index].cov = cov;
            config.arc_array[index].next = config.edge_array[contigId].arc;
            config.edge_array[contigId].arc = &config.arc_array[index];
            index ++ ;
        }
    }
    delete in1;
    assert(index == config.arcNum +1 );
}
void loadCluster(GlobalConfig & config)
{
    std::string line;
    unsigned int contigId;
    unsigned int to;
    float cov;
    config.connectionNum= 0;
    // Counting arcs
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(config.updateEdge);
    while(!std::getline(*in,line).eof())
    {
        std::istringstream ist(line);
        ist>>contigId;
        config.keys.insert(contigId);
        while(! ist.eof() )
        {
            ist>>to>>cov;
            config.connectionNum++ ;
            config.connections[contigId][to] = cov;
            config.connections[to][contigId] = cov ;
            config.keys.insert(to);
        }
    }
    delete in ;
    for( const auto & i : config.keys)
        config.edge_array[i].SetKey();
}
void buildConnection(GlobalConfig & config )
{
    for( const auto & i : config.keys)
    {
        std::stack<Edge> stack;
        std::map<unsigned int , Edge > history;
        std::map<unsigned int , std::vector<std::stack<Edge> > > paths;
        std::map<unsigned int , std::vector<std::stack<Edge> > > mids;
        config.edge_array[i].DepthSearch( config.edge_array , stack,
                history, paths , mids ,0 , config.connections.at(i) );

        
    }
}

}
}
