#include "biocommon/sam_bam/sam_parser.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <cassert>

int main()
{
    std::string line ;
    while ( ! std::getline( std::cin , line ).eof() )
    {
        BGIQD::SAM::LineParser parse(line) ;
        if( parse.IsVaid() && ! parse.IsHead() )
        {
            auto match_data = parse.ParseAsMatchData() ;
            std::string query_name = match_data.read_name ;
            std::string ref_name = match_data.ref_name;
            std::cout<<query_name<<'\t'
                <<ref_name<<'\t';
            if( match_data.ref_name != "*" )
            {
                std::cout<<"0\t0\t0\t0\t0\t0\t0\n";
            }
            else
            {
                auto & info = match_data ;
                int total_read_len  = info.read_len ;
                std::cout<<total_read_len<<'\t';
                std::cout<<info.first_match_position<<'\t';
                int total_result_len = info.total_result_len() ;
                int idy = info.md_data.total_same ;
                int missmatch = ( info.total_match_len() - info.md_data.total_same ) ;
                int in_f =  info.total_in_len() ;
                int del_f = info.total_del_len() ;

                assert( idy + missmatch + in_f + del_f == total_result_len );
                int m_base =     (info.total_match_len()) ;
                int in_base =    (info.total_in_len()) ;
                int clip_base =  (info.total_clip_len()) ;
                assert( m_base + in_base + clip_base   == total_read_len );

                std::cout<<idy<<'\t'
                    <<m_base<<'\t'
                    <<in_f<<'\t'
                    <<del_f<<'\t'
                    <<clip_base<<'\n';
            }
        }
        else
        {
        }
    }
    return 0 ;
}
