#include "common/args/argsparser.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "stLFR/CBB.h"
#include <map>

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(int , contig_num , "the contig id ");
    END_PARSE_ARGS

    unsigned int target_contig_num = contig_num.to_int() ;
    int bin_max = 0 ;
    std::string line ;

    std::map< int , BGIQD::stLFR::BinRelation > bin_cache ;


    while( ! std::getline( std::cin , line ).eof() )
    {
        BGIQD::stLFR::BinRelation tmp ;
        tmp.InitFromString(line);
        if( tmp.contigId == target_contig_num )
        {
            bin_cache[tmp.binId] = tmp ;
            if( tmp.binId > bin_max )
                bin_max = tmp.binId ;
        }
    }

    float * array = (float*)malloc(sizeof(float) * bin_max * bin_max ) ;
    for( int i = 0 ; i< bin_max ; i ++ )
    {
        for( int j = 0 ; j < bin_max ; j++ )
        {
            if( i == j )
                array[i*bin_max+j] = 1.0f ;
            else 
            {
                array[i*bin_max+j] = 0.0f ;
            }
        }

    }
    for( const auto & pair1 : bin_cache )
    {
        int i = pair1.first -1 ;
        for( const auto & info : pair1.second.sims )
        {
            int j = info.second.binId -1 ;
            array[i*bin_max+j] = info.second.simularity  ;
        }
    }

    for( int i = 0 ; i< bin_max ; i ++ )
    {
        for( int j = 0 ; j < bin_max ; j++ )
        {
            std::cout<< array[i*bin_max+j] ;
            if( j < bin_max -1 )
                std::cout<<'\t';
            else
                std::cout<<'\n';
        }
    }
    delete array ;
    return 0 ;
}

