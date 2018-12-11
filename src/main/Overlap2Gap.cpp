#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"
#include "common/freq/freq.h"

#include "soap2/fileName.h"

#include "stLFR/ScaffInfo.h"
#include "stLFR/ContigOverlap.h"

#include "biocommon/paf/PAF.h"

#include <map>
#include <set>
#include <vector>

struct AppConfig
{

    BGIQD::LOG::logger loger;

    std::string overlaps;

    std::map<std::string , std::vector<BGIQD::stLFR::OverlapInfo> >  overlap_caches ;

    BGIQD::FREQ::Freq<int> cache_freq ;

    BGIQD::FREQ::Freq<int> used_freq ;

    void Init()
    {
        BGIQD::LOG::logfilter::singleton().get("ScaffInfo2Seqs",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadOverlaps()
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(overlaps);
        if ( in == NULL )
            FATAL( " failed to open overlap file to read !!! ");
        std::string line ;
        while(! std::getline(*in , line).eof() )
        {
            BGIQD::PAF::PAF_Item tmp ;
            tmp.InitFromString(line);
            BGIQD::stLFR::OverlapInfo ov;
            ov.InitFromPAF(tmp);
            overlap_caches[ov.KeyT()].push_back(ov);
        }
        delete in;
    }

    void PrintFreq()
    {
        for( const auto & pair : overlap_caches)
        {
            cache_freq.Touch(int(pair.second.size()));
        }
        std::cerr<<"Cache freq\n"<<cache_freq.ToString();
        std::cerr<<"Used freq\n"<<used_freq.ToString();
    }

    void FilterAllScaffItems()
    {
        std::string line ;
        BGIQD::stLFR::ContigDetail prev ;
        prev.scaff_id = -1 ;
        while( ! std::getline(std::cin , line).eof() )
        {
            if( line[0] == '>' )
            {
                if( prev.scaff_id != -1 )
                {
                    std::cout<<prev.ToString()<<'\n';
                }
                std::cout<<line<<'\n';
                continue ;
            }
            BGIQD::stLFR::ContigDetail tmp ;
            tmp.InitFromString(line);
            if( tmp.scaff_id == prev.scaff_id )
            {
                BGIQD::stLFR::OverlapInfo ov;
                ov.InitFromScaffItem(prev,tmp);
                auto itr = overlap_caches.find( ov.KeyT() );
                if( itr == overlap_caches.end() )
                {
                    used_freq.Touch(0);
                }
                else 
                {
                    used_freq.Touch(int(itr->second.size()));
                    auto overlap = itr->second[0] ;
                    prev.gap_size = - overlap.overlap_len ;
                }
                std::cout<<prev.ToString()<<'\n';
            }
            else
            {
                ;
            }
            prev = tmp ;
        }
        if( prev.scaff_id != -1 )
        {
            std::cout<<prev.ToString()<<'\n';
        }
    }

} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, overlaps,"the overlaps.paf file ");
    END_PARSE_ARGS;

    config.overlaps = overlaps.to_string() ;
    config.Init( );
    BGIQD::LOG::timer t(config.loger,"Overlap2Gap");

    config.LoadOverlaps();
    config.FilterAllScaffItems();
    config.PrintFreq();
    return 0;
}
