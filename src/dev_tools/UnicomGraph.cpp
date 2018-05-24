#include "contig_barcode.h"
#include "common/args/argsparser.h"
#include <iostream>
#include "algorithm/disjoin_set/disjoin_set.h"
using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main()
{
    std::string line ;
    BGIQD::Algorithm::DisJoin_Set<int> dj_set;
    std::map< int ,std::map<int,int>  > nodeContig;
    while(!std::getline(std::cin,line).eof())
    {
        int fromV , endV , contigId , contigLen;
        sscanf(line.c_str() ,"\tV%d -> V%d[label =\"%d(%d)\"];\n",&fromV , &endV , &contigId, &contigLen);
        dj_set.AddConnect(fromV,endV);
        nodeContig[fromV][endV] = contigId;
    }

    for( const auto & i : nodeContig )
    {
        for( const auto & j : i.second){
            std::cout<<j.second<<"\t"<<dj_set.GetGroup(i.first)<<std::endl;
        }
    }

    return 0;
}
