#include "contig_barcode.h"
#include "argsparser.h"
#include "file_reader.h"
#include "stringtools.h"
#include "sam_parser.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <set>

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

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
void print_show(const std::string & line ,contigPos   & data , lines & l )
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
    updateMap1(l,contig,c1,f1);
    updateMap1(l,contig,c2,f2);
}

void loadGraph(lines & l , contigPos & c)
{
    BGIQD::LOG::timer t(BGIQD::log1,"loadGraph");
    std::string line ;
    while( !std::getline(std::cin, line).eof() )
    {
        print_show(line,c,l);
    }
}

void drawGraph( lines & l , contigPos & c)
{
    BGIQD::LOG::timer t(BGIQD::log1,"drawGraph");
    std::cout<<"graph{"<<std::endl;
    for(const auto & i : l)
    {
        std::cout<<i.first.first<<" -- "<<i.first.second
            <<"[label=\""<<i.second<<"\"]"
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
    loadGraph(l,p);
    drawGraph(l,p);
    return 0;
}
