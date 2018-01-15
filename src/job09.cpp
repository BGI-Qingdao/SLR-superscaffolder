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
typedef std::vector< std::tuple< int , int , float , int > > cluter_show;

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

typedef std::map<int , std::vector< std::tuple<int ,int> > > contigRef;
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
        data[read].push_back(std::make_tuple(d0.first_match_position,end_pos));
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

void print_show( int seed, const cluter_show & data)
{
    std::cout<<seed<<":";
    for( const auto & i : data )
    {
        std::cout<<"\t"<<std::get<0>(i)<<"-"<<std::get<1>(i)<<"("<<std::get<3>(i)<<":"<<std::get<2>(i)<<")";
    }
    std::cout<<std::endl;
    return ;
}

void print_contigPos(const contigRef & data)
{
    std::vector<std::tuple<int , int , int >> ret ;
    for( const auto  & a : data )
    {
        for( const auto & b : a.second)
        {
            ret.push_back(std::make_tuple( std::get<0>(b) , std::get<1>(b) , a.first ));
        }
    }
    std::sort( ret.begin() , ret.end() );
    for( const auto & c : ret )
    {
        if ( std::get<1>(c) > std::get<0>(c) )
        std::cout<<std::get<0>(c)<<"\t"
                << std::get<1>(c)<<"\t"
                << std::get<1>(c) - std::get<0>(c)<<"\t"
                << std::get<2>(c)<<std::endl;
    }
}

int main(int argc , char ** argv)
{
    initLog("JOB09");

    START_PARSE_ARGS
    DEFINE_ARG(std::string , refBarcode , 'i');
    DEFINE_ARG(std::string , refContig, 'c');
    DEFINE_ARG(std::string , output, 'o');
    DEFINE_ARG(bool , ponly, 'p');
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
        print_show(i.first, get_show( i.second, r ));
    }
    return 0;
}
