#include "algorithm/collection/collection.h"
#include "algorithm/interval/Interval.h"

#include "common/args/argsparser.h"
#include "common/freq/freq.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"

#include <set>
#include <map>
#include <vector>

struct GlobalConfig
{
    typedef BGIQD::Collection::Collection<unsigned int> SeedVec;
    typedef BGIQD::INTERVAL::Interval<int> gapLen;
    struct GapInfo
    {
        unsigned int prev;
        unsigned int next ;
        BGIQD::Collection::Collection<unsigned int>  seedsInRef;
        BGIQD::Collection::Collection<unsigned int>  seedsCluster;
    };

    std::string seedLinear ;
    std::string seedCluster;

    std::map<unsigned int , std::set<int> > seedPos;
    std::map<unsigned int , int > seedLens;
    std::map<int,unsigned int> pos2Contig;
    std::vector<GapInfo> gaps;
    BGIQD::FREQ::Freq<int>  freq;
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
            pos2Contig[line_index] = std::stoul(items[0]);
            line_index ++ ;
        }
        delete sf ;
    }

    void LoadSeedCluster()
    {
        auto sf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(seedCluster);

        std::string line;

        while(!std::getline(*sf,line).eof())
        {
            auto items = BGIQD::STRING::split(line,"\t");
            GapInfo info ;
            info.prev = std::stoul(items[0]);
            info.next = std::stoul(items[1]);

            for( int i = 0 ; i < (int) items.size(); i++ )
            {
                info.seedsCluster.IncreaseElement(std::stoul(items[i]));
            }

            gapLen min(0,0);

            for( auto x : seedPos[info.prev])
            {
                gapLen len ;
                len.min = x ;
                for( auto y : seedPos[info.next] )
                {
                    len.max = y ;
                    if( min.Len() == 0 || std::abs(len.Len()) < std::abs(min.Len()) )
                    {
                        min = len ;
                    }
                }
            }
            if( std::abs( min.Len()) <= 1 )
            {
                freq.Touch(101);
                continue ;
            }
            int start = std::min(min.max, min.min ) +1 ;
            int end= std::max(min.max, min.min ) ;
            for( int i = start ; i < end ; i ++ )
            {
                info.seedsInRef.IncreaseElement(pos2Contig[i]);
            }
            auto both = SeedVec::Intersection(info.seedsCluster,info.seedsInRef);
            freq.Touch(both.keysize()*100/info.seedsInRef.keysize());
        }
        std::cerr<<freq.ToString()<<std::endl;
        delete sf ;
    }

} config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS 
    DEFINE_ARG_REQUIRED(std::string , seedLinear , " the seed linear file ");
    DEFINE_ARG_REQUIRED(std::string , gap_cluster, " the gap cluster file ");
    END_PARSE_ARGS

    config.seedLinear = seedLinear.to_string();
    config.seedCluster = gap_cluster.to_string() ;
    config.LoadSeedLinear();
    config.LoadSeedCluster();
    return 0;
}
