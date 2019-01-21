#include "biocommon/paf/PAF.h"
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
        BGIQD::PAF::PAF_Item item ;
        item.InitFromString(line) ;
        std::string query_name = item.query_name;
        std::string ref_name = item.target_name;
        std::cout<<query_name<<'\t'
            <<ref_name<<'\t';
        if( ref_name == "*" )
        {
            std::cout<<"0\t0\t0\t0\t0\t0\t0\n";
        }
        else
        {
            auto & info = item.details;
            int total_read_len  = item.query_len ;
            std::cout<<total_read_len<<'\t';
            std::cout<<item.target_start<<'\t';
            int total_result_len = info.total_result_len() ;
            int idy = item.md_data.total_same ;
            int missmatch = ( info.total_match_len() - idy) ;
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
    return 0 ;
}
