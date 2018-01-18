#include "contig_barcode.h"
#include "argsparser.h"
#include "file_reader.h"
#include "stringtools.h"
#include "sam_parser.h"
#include <algorithm>
#include <iostream>
#include <sstream>

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

typedef std::map< int , std::map< int , float > >  cluters;
typedef std::vector< std::tuple< int , int , float , int > > cluter_show;

void updateMap1(std::map<int,float> & map , int key , float value )
{
    auto itr = map.find(key) ;
    if( itr == map.end() ||  itr->second < value ) 
        map[key] = value ;
}

void incrMap1(std::map<int,int> & map , int key  )
{
    auto itr = map.find(key) ;
    if( itr == map.end() ) 
        map[key] = 1 ;
    else
        itr->second ++ ;
}
void print_show(const std::string & line , std::map<int , int > & data )
{
    std::istringstream ism(line);
    int contig , gap ;
    ism>>contig>>gap;
    int pend = -1 ;
    while(!ism.eof())
    {
        int start , end , c , f ; 
        char cc;
           // 101   -   104  (  10  :   0.44 )
        ism>>start>>cc>>end>>cc>>c>>cc>>f>>c;
        if(pend != -1 )
        {
            incrMap1(data, start-pend );
            pend = end ;
        }
    }
}
int main()
{
    std::map< int , int > gaps;
    initLog("JOB10");
    {
        std::string line ;
        while( !std::getline(std::cin, line).eof() )
        {
            print_show(line,gaps);
        }
    }
    for(const auto & i : gaps )
    {
        std::cout<<i.first<<"\t"<<i.second<<std::endl;
    }

    return 0;
}
