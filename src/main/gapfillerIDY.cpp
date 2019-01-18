#include "biocommon/sam_bam/sam_parser.h"
#include <string>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <set>

std::map<long , std::vector<BGIQD::SAM::MatchData> >   cache_data;

std::set<long> querys;
std::set<long> matched;

int main()
{
    std::string line ;
    while ( ! std::getline( std::cin , line ).eof() )
    {
        BGIQD::SAM::LineParser parse(line) ;
        if( parse.IsVaid() && ! parse.IsHead() )
        {
            auto match_data = parse.ParseAsMatchData() ;
            long query_num = std::stoi(match_data.read_name) ;
            querys.insert(query_num);
            if( match_data.ref_name != "*" )
                matched.insert(query_num);
            if( match_data.ref_name != match_data.read_name )
                continue ;
            cache_data[query_num].push_back(match_data);
        }
    }

    auto less = []( const BGIQD::SAM::MatchData & r1 ,
            const BGIQD::SAM::MatchData & r2 ) -> bool
    {
        return r1.total_result_len() < r2.total_result_len() ;
    };

    for( auto & pair : cache_data )
    {
        if( pair.second.size() > 1 )
            std::sort( pair.second.rbegin() , pair.second.rend() , less );
    }
    std::cerr<<"total query   "<<querys.size()<<std::endl;
    std::cerr<<"total matched "<<matched.size()<<std::endl;
    std::cerr<<"total self matched "<<cache_data.size()<<std::endl;

    long idy = 0 ;
    long missmatch = 0 ;
    long in_f = 0;
    long del_f = 0;

    long  m_base = 0 ;
    long  in_base = 0 ;
    long  clip_base = 0 ;
    long total_result_len = 0 ;
    long total_read_len = 0;
    for( const auto & pair : cache_data )
    {
        const auto & info = pair.second.at(0);
        idy+= info.md_data.total_same ;
        total_result_len += info.total_result_len() ;
        total_read_len  += info.read_len ;
        missmatch += ( info.total_match_len() - info.md_data.total_same ) ;
        in_f =  info.total_in_len() ;
        del_f = info.total_del_len() ;

        m_base +=     (pair.second[0].total_match_len()) ;
        in_base +=    (pair.second[0].total_indel_len()) ;
        clip_base +=  (pair.second[0].total_clip_len()) ;
    }
    int total_size = cache_data.size() ;
    std::cerr<<"mean idy     " <<float(idy) / float(total_result_len) / float(total_size)<<std::endl;
    std::cerr<<"mean missmatch     " << (float)missmatch / float(total_result_len) / float(total_size) <<std::endl;
    std::cerr<<"mean insert " << (float)in_f / float(total_result_len) / float(total_size) <<std::endl;
    std::cerr<<"mean deletion " << (float)del_f  / float(total_result_len) / float(total_size) <<std::endl;
    std::cerr<<"________________"<<std::endl;
    std::cerr<<"mean match fraction " << (float) m_base / float(total_read_len) / float(total_size) <<std::endl;
    std::cerr<<"mean insert fraction " << (float)in_base  / float(total_read_len) / float(total_size) <<std::endl;
    std::cerr<<"mean clip fraction " << (float) clip_base / float(total_read_len) / float(total_size) <<std::endl;

    return 0 ;
}
