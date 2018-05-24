#include "contig_barcode.h"
#include "argsparser.h"
#include "file_reader.h"
#include "stringtools.h"
#include "sam_parser.h"
#include <algorithm>
#include <iostream>
using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

typedef std::map< int , std::map< int , float > >  cluters;

void updateMap1(std::map<int,float> & map , int key , float value )
{
    auto itr = map.find(key) ;
    if( itr == map.end() ||  itr->second < value ) 
        map[key] = value ;
}

void loadCluterData( cluters & data)
{
    //auto  in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
    std::string line ;

    while(! getline(std::cin,line).eof() )
    {
        auto dd = BGIQD::STRING::split(line,"\t");
        int key = 0 ;
        for(const auto i:dd)
        {
            auto o = BGIQD::STRING::split(i,":");
            if(o[1] == "1")
            {    key = std::stoi(o[0]); break;}
        }
        for(const auto i:dd)
        {
            auto o = BGIQD::STRING::split(i,":");
            updateMap1(data[key],std::stoi(o[0]),std::stof(o[1]));
        }
    }
    //delete in ;
}



int main()
{
    initLog("MergeClusterResult");

    cluters c;
    loadCluterData( c);
    for( const auto i : c )
    {
        int key = i.first;
        if( i.second.empty() )
            continue ;
        std::cout<<i.first;
        for( const auto j : i.second)
        {
            if( j.first != key )
            {
                std::cout<<"\t"<<j.first<<"\t"<<j.second;
            }
        }
        std::cout<<std::endl;
    }
    return 0;
}
