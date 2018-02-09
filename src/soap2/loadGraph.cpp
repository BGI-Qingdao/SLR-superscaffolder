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
    sscanf(line.c_str(),"EDGEs %u",&config.contigTotalNum);
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
        if( length > config.K )
            config.edge_array[index].length = length;
        else
            config.edge_array[index].length = 0 ;
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
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(config.updateEdge);
    // load connection 
    while(!std::getline(*in,line).eof())
    {
        std::istringstream ist(line);
        ist>>contigId;
        config.keys.insert(contigId);
        while(! ist.eof() )
        {
            ist>>to>>cov;
            config.connections[contigId][to] = cov;
            config.connections[to][contigId] = cov ;
            config.keys.insert(to);
        }
    }
    delete in ;

    // init keys
    config.clusterNum = config.keys.size();
    config.key_array =static_cast<KeyEdge*>( calloc( sizeof(KeyEdge*) , config.clusterNum +1));
    unsigned int index = 1;
    for( const auto & i : config.keys)
    {
        config.edge_array[i].SetKey();
        config.key_array[index].edge_id = i;
        config.key_array[index].id = index;
        config.key_map[i] = index ;
        index++;
    }
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
                history, paths , mids ,config.edge_array[i].length , config.connections.at(i) );

        for(const auto & j : paths)
        {
            config.key_array[i].to.insert(j.first);
            config.key_array[j.first].from.insert(i);
        }
    }
}

void LinearConnection(GlobalConfig &config)
{
    for( unsigned int i = 1 ; i< config.clusterNum +1 ; i++ )
    {
        KeyEdge & curr =  config.key_array[i];
        if( curr.IsMarked() 
            ||curr.IsLinear() 
            || curr.IsTipTo()
            )
            continue;
        if( curr.IsSingle() )
        {
            std::vector<unsigned int > a;
            a.push_back(curr.edge_id) ;
            config.contigs.push_back(a);
            curr.Mark();
        }
        else
        {
            for(auto next : curr.to )
            {
                std::vector<unsigned int > path;
                path.push_back(next) ;
                unsigned int next_k = config.key_map[next];
                while( config.key_array[next_k].IsLinear() )
                {
                    unsigned int next_i = *config.key_array[next_k].to.begin();
                    next_k = config.key_map[next_i];
                    path.push_back(next_i) ;
                }
                unsigned int next_e = *config.key_array[next_k].to.begin();
                path.push_back(next_e) ;
                config.contigs.push_back(path);
            }
        }
    }
}


}
}
