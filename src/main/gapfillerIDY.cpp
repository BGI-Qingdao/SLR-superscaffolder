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
    float idy_total = 0.0f ;
    float err_base = 0.0f ;
    for( const auto & pair : cache_data )
    {
        idy_total += pair.second[0].md_data.IDY(pair.second[0].total_match_len());
        err_base +=  (float)(pair.second[0].total_match_len()) / (float)(pair.second[0].read_len) ;
    }
    idy_total = idy_total / cache_data.size() ;
    err_base = err_base / cache_data.size() ;
    std::cerr<<"mean idy     " <<idy_total<<std::endl;
    std::cerr<<"mean match fraction " <<err_base<<std::endl;

    return 0 ;
}
