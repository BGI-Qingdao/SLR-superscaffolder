#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"
#include "common/freq/freq.h"

#include "algorithm/interval/Interval.h"

#include <iostream>
#include <set>

struct AppConfig
{
    typedef BGIQD::INTERVAL::Interval<int,BGIQD::INTERVAL::IntervalType::Left_Close_Right_Close> ContigArea;

    std::string seedLinear ;
    std::string gap;

    std::map<unsigned int , std::vector<ContigArea> > seedArea;
    std::map<unsigned int , int > seedLens;

    void LoadSeedLinear()
    {
        auto sf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(seedLinear);
        std::string line;
        int line_index = 1 ;
        while(!std::getline(*sf,line).eof())
        {
            auto items = BGIQD::STRING::split(line,"\t");
            seedArea[std::stoul(items[0])].push_back(ContigArea(std::stoi(items[2]),std::stoi(items[3])));
            seedArea[std::stoul(items[0])+1].push_back(ContigArea(std::stoi(items[2]),std::stoi(items[3])));
            line_index ++ ;
        }
        delete sf ;
    }

    void ParseGaps()
    {
        auto sf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(gap);
        std::string line;
        while(!std::getline(*sf,line).eof())
        {
            auto items = BGIQD::STRING::split(line,"\t");
            unsigned int p = std::stoi(items[0]);
            unsigned int e = std::stoi(items[1]);

            if(seedArea[p].size() != 1 || seedArea[e].size() != 1 )
            {
                std::cout<<p<<'\t'<<e<<'\t'<<'*'<<'\n';
            }
            else 
            {
                auto & p_area = seedArea[p][0];
                auto & e_area = seedArea[e][0];
                if( p_area.min < e_area.min )
                {
                    std::cout<<p<<'\t'<<e<<'\t'<<e_area.min - p_area.max<<'\n';
                }
                else
                {
                    std::cout<<p<<'\t'<<e<<'\t'<<p_area.min - e_area.max<<'\n';
                }
            }
        }
    }
}config;


int main(int argc , char **argv)
{
    START_PARSE_ARGS 
    DEFINE_ARG_REQUIRED(std::string , seedLinear , " the seed linear file ");
    DEFINE_ARG_REQUIRED(std::string , gap, " the trunk linear file ");
    END_PARSE_ARGS

    config.seedLinear = seedLinear.to_string() ;
    config.gap = gap.to_string() ;

    config.LoadSeedLinear();
    config.ParseGaps();
    return 0;
}
