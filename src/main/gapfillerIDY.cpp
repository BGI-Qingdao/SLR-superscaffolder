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
    float idy = 0.0f ;
    float missmatch = 0.0f ;
    float in_f ;
    float del_f ;
    
    float m_base = 0.0f ;
    float in_base = 0.0f ;
    float clip_base = 0.0f ;
    for( const auto & pair : cache_data )
    {
        const auto & info = pair.second.at(0);
        idy+= info.md_data.IDY(info.total_match_len());
        missmatch +=  float(info.total_match_len() - info.md_data.total_same ) /(float) info.total_result_len();
        in_f = float( info.total_in_len() ) / (float) info.total_result_len();
        del_f = float( info.total_del_len() ) /(float) info.total_result_len();

        m_base +=  (float)(pair.second[0].total_match_len()) / (float)(pair.second[0].read_len) ;
        in_base +=  (float)(pair.second[0].total_indel_len()) / (float)(pair.second[0].read_len) ;
        clip_base +=  (float)(pair.second[0].total_clip_len()) / (float)(pair.second[0].read_len) ;
    }
    idy= idy/ cache_data.size() ;
    clip_base= clip_base/ cache_data.size() ;
    in_base= in_base/ cache_data.size() ;
    missmatch= missmatch/ cache_data.size() ;
    in_f= in_f/ cache_data.size() ;
    del_f= del_f/ cache_data.size() ;
    m_base = m_base / cache_data.size() ;

    std::cerr<<"mean idy     " <<idy<<std::endl;
    std::cerr<<"mean missmatch     " <<missmatch<<std::endl;
    std::cerr<<"mean insert " <<in_f<<std::endl;
    std::cerr<<"mean deletion " <<del_f<<std::endl;
    std::cerr<<"________________"<<std::endl;
    std::cerr<<"mean match fraction " <<m_base<<std::endl;
    std::cerr<<"mean insert fraction " <<in_base<<std::endl;
    std::cerr<<"mean clip fraction " <<clip_base<<std::endl;

    return 0 ;
}
