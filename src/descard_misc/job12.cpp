#include "contig_barcode.h"
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "biocommon/sam_bam/sam_parser.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <set>
#include <algorithm/disjoin_set/disjoin_set.h>

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;
using namespace BGIQD::Algorithm;

typedef std::map< int , std::map< int , int> >  contigPos;
typedef std::map<std::pair<int , int > ,float> lines;

void updateMap1(lines&map1, int key1 , int key2, float v)
{
    int k1 = key1 < key2 ? key1 : key2 ;
    int k2 = key1 < key2 ? key2 : key1 ;
    auto itr = map1.find(std::make_pair(k1,k2)) ;
    if( itr == map1.end() ) 
        map1[std::make_pair(k1,k2)] =v;
}

void incrMap1(std::map<int,int> & map , int key  )
{
    auto itr = map.find(key) ;
    if( itr == map.end() ) 
        map[key] = 1 ;
    else
        itr->second ++ ;
}
void print_show(const std::string & line ,contigPos   & data , lines & l ,DisJoin_Set<int> &d)
{
    std::istringstream ism(line);
    int contig , gap ;
    ism>>contig>>gap;
    int c1 , c2 ;
    float f1 =0.0f , f2 = 0.0f;
    while(!ism.eof())
    {
        int start , end , c ;
        float f ; 
        char cc;
           // 101   -   104  (  10  :   0.44 )
        ism>>start>>cc>>end>>cc>>c>>cc>>f>>cc;
        data[c][start]=end;
        if ( f == 1.0f )
            continue ;
        if(  f < 0.14f )
            continue;
        if( f > f1 )
        {
            f1 = f ;
            c1 = c;
        }
        else if ( f > f2 )
        {
            f2 = f ;
            c2 = c ;
        }
    }
    if ( f1 > 0.14f )
    {
    updateMap1(l,contig,c1,f1);
    //updateMap1(l,contig,c2,f2);
    d.AddConnect(c1,contig);
    }
}

void loadGraph(lines & l , contigPos & c , DisJoin_Set<int> & d)
{
    BGIQD::LOG::timer t(BGIQD::log1,"loadGraph");
    std::string line ;
    while( !std::getline(std::cin, line).eof() )
    {
        print_show(line,c,l,d);
    }
}

void drawGraph( lines & l , contigPos & c , DisJoin_Set<int> & d)
{
    BGIQD::LOG::timer t(BGIQD::log1,"drawGraph");
    std::cout<<"graph{"<<std::endl;
    for(const auto & i : l)
    {
        std::cout<<i.first.first<<" -- "<<i.first.second
            <<"[label=\""<<i.second<<"\"("
            <<d.GetGroup(i.first.first)<<")]"
            <<std::endl;
    }
    std::cout<<"}"<<std::endl;
}

int main()
{
    std::map< int , int > gaps;
    initLog("JOB10");
    contigPos p;
    lines l;
    DisJoin_Set<int> d;
    loadGraph(l,p,d);
    drawGraph(l,p,d);
    return 0;
}
