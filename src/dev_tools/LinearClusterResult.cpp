#include "contig_barcode.h"
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "biocommon/sam_bam/sam_parser.h"
#include <algorithm>
#include <iostream>
using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

typedef std::map< int , std::map< int , float > >  cluters;
typedef std::vector< std::tuple< int , int , float , int > > cluter_show;

void incrMap2(std::map<int,int> & map , int key )
{
    auto itr = map.find(key) ;
    if( itr == map.end()  ) 
        map[key] = 1 ;
    else
        itr->second ++ ;
}
void updateMap1(std::map<int,float> & map , int key , float value )
{
    auto itr = map.find(key) ;
    if( itr == map.end() ||  itr->second < value ) 
        map[key] = value ;
}

void loadCluterData(const std::string & file , cluters & data)
{
    auto  in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
    std::string line ;

    while(! getline((*in),line).eof() )
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
    delete in ;
}

typedef std::map<int , std::vector< std::tuple<int ,int,bool> > > contigRef;
void loadContigRef( const std::string & file , contigRef & data )
{
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
    while( ! in->eof() )
    {
        std::string line;
        std::getline(*in,line);
        if( in->eof() )
            break;
        BGIQD::SAM::LineParser p(line);
        if(p.IsHead())
            continue;
        auto d0 = p.ParseAsMatchData();
        int read = std::stoi(d0.read_name);
        int end_pos = 0;
        for( size_t i = 0 ; i< d0.detail.infos.size() ; i++ )
        {
            auto info = d0.detail.infos[i];
            if( info.type != BGIQD::SAM::M )
            {
                continue;
            }
            end_pos = info.end_position_on_ref ;
        }
        data[read].push_back(std::make_tuple(d0.first_match_position,end_pos, d0.IsReverseComplete()));
    }
    delete in ;
}

cluter_show get_show( const std::map<int , float>d1 , const contigRef & data)
{
    cluter_show ret ;
    for ( const auto & i : d1 )
    {
        auto j = data.at(i.first);
        for( const auto jj : j)
        {
            ret.push_back(std::make_tuple( std::get<0>(jj) , std::get<1>(jj) , i.second , i.first));
        }
    }
    std::sort(ret.begin() , ret.end());
    return ret;
}

cluter_show filter_repeat_in_region(const cluter_show & data ,bool repeat, int from, int to)
{
    std::map<int,int> freq;
    std::set<int> mayerr;
    int key = -1 ;
    cluter_show newer;
    for ( size_t i = 0 ; i < data.size() ; i++ )
    {
        incrMap2(freq,std::get<2>(data[i]));
        if( i >=1 )
        {
            if( std::get<0>(data[i]) - std::get<1>(data[i-1]) > 1000000 )
            {
                mayerr.insert(i);
                mayerr.insert(i-1);
            }
        }
        if( from >-1 && to > 0 )
        {
            if( std::get<2>(data[i] ) - 1.0f >=-0.00001 )
            {
                key = std::get<3>(data[i]);
                if(! ( (std::get<0>(data[i])>=from && std::get<0>(data[i])<=to) || (std::get<1>(data[i])>=from && std::get<1>(data[i])<=to) ))
                    return newer;
            }
        }
    }
    if( repeat )
        return data;
    bool keyfound = false ;
    for( size_t i = 0 ; i < data.size() ; i++ )
    {
        if ( mayerr.find(i) != mayerr.end()  && freq.at(std::get<2>(data[i])) > 1 )
        {
            if( keyfound )
                break;
            else
                newer.clear();
            continue;
        }
        if ( std::get<3>(data[i]) == key )
            keyfound  = true;  
        newer.push_back(data[i]);
    }
    return newer;
}

void print_show( int seed, const cluter_show & data , bool s )
{
    if ( data.size() < 2 )
        return ;
    int first = std::get<0>(*data.begin());
    int end = std::get<0>(*data.rbegin());
    std::cout<<seed<<"\t"<<end-first;
    if(s)
        std::cout<<"\t"<<first<<"\t"<<end;
    for( const auto & i : data )
    {
        std::cout<<"\t"<<std::get<0>(i)<<"-"<<std::get<1>(i)<<"("<<std::get<3>(i)<<":"<<std::get<2>(i)<<")";
    }
    std::cout<<std::endl;
    return ;
}

void print_contigPos(const contigRef & data)
{
    std::vector<std::tuple<int , int , int , bool >> ret ;
    for( const auto  & a : data )
    {
        for( const auto & b : a.second)
        {
            ret.push_back(std::make_tuple( std::get<0>(b) , std::get<1>(b) , a.first, std::get<2>(b)));
        }
    }

    std::sort( ret.begin() , ret.end() );

    for( const auto & c : ret )
    {
        if ( std::get<1>(c) > std::get<0>(c) )
        std::cout<<std::get<0>(c)<<"\t"
                << std::get<1>(c)<<"\t"
                << std::get<1>(c) - std::get<0>(c)<<"\t"
                << std::get<2>(c) <<'\t'
                << std::get<3>(c) <<std::endl;
    }

}

int main(int argc , char ** argv)
{
    initLog("LinearClusterResult");

    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , refBarcode , "the clusters result");
    DEFINE_ARG_REQUIRED(std::string , refContig,"sam file . map contig to ref");
    DEFINE_ARG_OPTIONAL(bool , ponly,  "If this flag setted , it will only print the postion of each contig.","false");
    DEFINE_ARG_OPTIONAL(bool , repeat,  "If this flag setted , it will not delete repeat contig from cluster.", "false");
    DEFINE_ARG_OPTIONAL(int , from , "from . default[0]","0");
    DEFINE_ARG_OPTIONAL(int , to , " to. default [0] , means all position is valid", "0");
    DEFINE_ARG_OPTIONAL(bool, show,"If this flag setted , it will print total length region in column 3 , 4 .","false");
    END_PARSE_ARGS

    cluters c;
    contigRef r;
    loadContigRef(refContig.to_string() , r);
    if( ponly.to_bool() )
    {
        print_contigPos(r);
        return 0;
    }
    loadCluterData(refBarcode.to_string(), c);
    for( const auto i : c )
    {
        print_show(i.first, filter_repeat_in_region( get_show( i.second, r ),repeat.to_bool(),from.to_int() , to.to_int()),show.to_bool());
    }
    return 0;
}
