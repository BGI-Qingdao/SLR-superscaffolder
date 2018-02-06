#include "contig_barcode.h"
#include "argsparser.h"
#include <iostream>

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main()
{
    initLog("ContigGraphType");
    std::string line ;
    std::map<int , int> inP;
    std::map<int , int> outP;
    std::set<int> nodes;
    std::map<int , std::tuple<int , int , int > > contigs;

    while(!std::getline(std::cin,line).eof())
    {
        int fromV , endV , contigId , contigLen;
        sscanf(line.c_str() ,"\tV%d -> V%d[label =\"%d(%d)\"];\n",&fromV , &endV , &contigId, &contigLen);
        if( inP.find( endV ) != inP.end() )
        {
            inP[endV]++;
        }
        else
        {
            inP[endV] = 1 ;
        }
        if( outP.find( fromV ) != outP.end() )
        {
            outP[fromV]++;
        }
        else
        {
            outP[fromV]=1;
        }
        nodes.insert(fromV);
        nodes.insert(endV);
        contigs[contigId] = std::make_tuple( fromV , endV , contigLen) ;
    }

    int zero = 0;
    int one = 0 ;
    int two = 0;
    std::map<int , int> zm;
    std::map<int , int> om;
    std::map<int , int> tm;

    for ( const auto & c : contigs ) 
    {
        int len ,fromV , endV;
        int index ;
        std::tie(fromV,endV,len) = c.second;
        //int fromV = std::get<0>(c.second) ;
        //int endV = std::get<1>(c.second);
        if( len <100)
            index = len/10 *10 ;
        else if ( len < 1000 )
            index = len / 100 * 100 ;
        else 
            index = len /1000 * 1000;
        if( /*inP.at( endV )==1 
        && outP.at( fromV ) ==1 
        && */inP.find(fromV) == inP.end() 
        && outP.find(endV) == outP.end() )
        {
            if( zm.find(index) == zm.end() )
                zm[index] = 0;
            zm[index] ++;
            std::cout<<c.first<<"\t"<<std::get<2>(c.second)<<"\t0"<<std::endl;
            zero ++ ;
        }
        else if (/*inP.at(endV) == 1  
        && outP.at( fromV ) ==1 
        && */( inP.find(fromV) == inP.end() 
        || outP.find(endV) == outP.end() ))
        {
            if( om.find(index) == om.end() )
                om[index] = 0;
            om[index] ++;
            std::cout<<c.first<<"\t"<<std::get<2>(c.second)<<"\t1"<<std::endl;
            one ++;
        }
        else
        {
            if( tm.find(index) == tm.end() )
                tm[index] = 0;
            tm[index] ++;
            std::cout<<c.first<<"\t"<<std::get<2>(c.second)<<"\t2"<<std::endl;
            two ++ ;
        }
    }
    std::cout<<"total nodes "<<nodes.size()<<std::endl;
    std::cout<<"from nodes "<<outP.size()<<std::endl;
    std::cout<<"end nodes "<<inP.size()<<std::endl;
    std::cout<<"contig "<<contigs.size()<<std::endl;
    std::cout<<" zero "<<zero<<" one "<<one<<" two "<<two<<std::endl;

    for( const auto& z : zm)
    {
        std::cout<<z.first<<"\t"<<z.second<<std::endl;
    }
    for( const auto& z : om)
    {
        std::cout<<z.first<<"\t"<<z.second<<std::endl;
    }
    for( const auto& z : tm)
    {
        std::cout<<z.first<<"\t"<<z.second<<std::endl;
    }
    return 0;
}
