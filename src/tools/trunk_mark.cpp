#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"
#include "common/freq/freq.h"

#include "algorithm/interval/Interval.h"

#include <iostream>
#include <set>

typedef BGIQD::INTERVAL::Interval<int,BGIQD::INTERVAL::IntervalType::Left_Close_Right_Close> ContigArea;
struct GlobalConfig 
{
    std::string seedLinear ;
    std::string trunkLinear;
    std::string outPrefix;
    int bin;
    std::map<unsigned int , std::set<int> > seedPos;
    std::map<unsigned int , std::set<ContigArea> > seedArea;
    std::map<unsigned int , int > seedLens;

    std::set<unsigned int> repeats;
    std::vector<ContigArea> scaffold_area;
    std::set<unsigned int> seeds ;
    void LoadSeedLinear()
    {
        auto sf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(seedLinear);
        std::string line;
        int line_index = 1 ;
        while(!std::getline(*sf,line).eof())
        {
            auto items = BGIQD::STRING::split(line,"\t");
            seedPos[std::stoul(items[0])].insert( line_index) ;
            seedLens[std::stoul(items[0])] = std::stoul(items[1]);
            seedArea[std::stoul(items[0])].insert(ContigArea(std::stoi(items[2]),std::stoi(items[3])));
            line_index ++ ;
        }
        delete sf ;
    }

    void ParseTrunkLinear()
    {
        auto tf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(trunkLinear);
        std::ostream * out = NULL ;
        int trunk_index = 1 ;

        std::string line ;
        ContigArea a_scaf(-1,-1);
        while(!std::getline(*tf,line).eof())
        {
            if(line[0] == '-')
            {
                if( out != NULL )
                    delete out ;
                out = NULL ;
                std::string of = outPrefix + std::to_string(trunk_index);
                trunk_index ++ ;
                out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(of);
                if( a_scaf.min != -1 && a_scaf.max != -1 )
                { 
                    if( a_scaf.min > a_scaf.max )
                        std::swap( a_scaf.min , a_scaf.max );
                    scaffold_area.push_back(a_scaf);
                }
                a_scaf = ContigArea(-1,-1);
                continue;
            }
            unsigned int contig = std::stoul(line);
            (*out)<<contig;
            if( seedPos[contig].size() > 1 )
            {/*
                if( a_scaf.min != -1 && a_scaf.max != -1 )
                {
                    if( a_scaf.min > a_scaf.max )
                        std::swap( a_scaf.min , a_scaf.max );
                    scaffold_area.push_back(a_scaf);
                }
                a_scaf = ContigArea(-1,-1);
                repeats.insert(contig);*/
            }
            else
            {
                seeds.insert( contig );
                if( a_scaf.min == -1 )
                {
                    a_scaf.min = seedArea[contig].begin()->min ;
                    a_scaf.max = seedArea[contig].begin()->max ;
                }
                else
                    a_scaf.max = seedArea[contig].begin()->min ;
            }

            for( auto x : seedPos[contig])
            {
                (*out)<<'\t'<<x;
            }
            (*out)<<'\n';
        }
        if( a_scaf.min != -1 && a_scaf.max != -1 )
        { 
            if( a_scaf.min > a_scaf.max )
                std::swap( a_scaf.min , a_scaf.max );
            scaffold_area.push_back(a_scaf);
        }
        delete tf ;
    }

    void Report()
    {

        std::cerr<<"____________________"<<std::endl;
        long total = 0;
        for(const auto & scaff : scaffold_area )
        {
            std::cerr<<scaff.min<<'\t'<<scaff.max<<'\t'<<scaff.Len()<<std::endl;
            total += scaff.Len() ;
        }

        std::cerr<<"____________________"<<std::endl;
        for(int i = 0 ; i <= 60000000 ; i+= bin)
        {
            ContigArea area(i , i+bin );
            int scaff_len = 0 ;
            int contig_len = 0;
            for( const auto & scaff : scaffold_area )
            {
                scaff_len += area.Overlap(scaff).Len();
            }

            for( const auto & s: seeds)
            {
                if( seedArea[s].size() == 1 )
                    contig_len += area.Overlap(*seedArea[s].begin()).Len();
            }
            std::cerr<<i<<'\t'<<scaff_len<<'\t'<<contig_len<<std::endl;
        }
        std::cerr<<"____________________"<<std::endl;
        long over = 0 ;
        for( int i = 0 ; i < (int)scaffold_area.size() ; i ++ )
        {
            for( int j = i + 1 ; j < (int)scaffold_area.size() ; j ++ )
            {
                auto overarea = scaffold_area[i].Overlap(scaffold_area[j]);
                over += overarea.Len() ;
            }
        }
        long contigLen = 0;
        for( const auto & s: seeds)
        {
            if( seedArea[s].size() == 1 )
            {
                contigLen += seedArea[s].begin()->Len();
            }
        }
        std::cerr<<"Total trunk len is "<<total<<std::endl;
        std::cerr<<"Total overlap len is "<<over<<std::endl;
        std::cerr<<"Total contig len is "<<contigLen<<std::endl;
    }

    void Init(const std::string & s , const std::string & t , const std::string & o)
    {
        seedLinear = s ;
        trunkLinear = t ;
        outPrefix = o ;
    }

}config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS 
    DEFINE_ARG_REQUIRED(std::string , seedLinear , " the seed linear file ");
    DEFINE_ARG_REQUIRED(std::string , trunkLinear , " the trunk linear file ");
    DEFINE_ARG_REQUIRED(std::string , outPrefix, " the out file prefix");
    DEFINE_ARG_OPTIONAL(int ,bin, " the bin size for calc freq","100000");
    END_PARSE_ARGS

    config.Init(seedLinear.to_string() , trunkLinear.to_string() , outPrefix.to_string() );
    config.bin = bin.to_int();
    config.LoadSeedLinear() ;
    config.ParseTrunkLinear();
    config.Report() ;
    return 0 ;
}
