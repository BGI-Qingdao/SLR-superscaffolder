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
    BGIQD::SOAP2::FileNames fNames;

    BGIQD::LOG::logger loger;

    std::string overlaps;

    std::map<std::string , std::vector<BGIQD::stLFR::OverlapInfo> >  overlap_caches ;

    BGIQD::FREQ::Freq<int> cache_freq ;

    BGIQD::FREQ::Freq<int> used_freq ;

    void Init( const std::string & prefix )
    {
        fNames.Init(prefix);
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

    void PrintCacheFreq()
    {
        for( const auto & pair : overlap_caches)
        {
            cache_freq.Touch(int(pair.second.size()));
        }
        std::cerr<<"Cache freq\n"<<cache_freq.ToString();
    }

} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.scaff_infos; Output xxx.scaff_infos_new");
        DEFINE_ARG_REQUIRED(std::string, overlaps,"the overlaps.paf file ");
    END_PARSE_ARGS;

    config.overlaps = overlaps.to_string() ;
    config.Init( prefix.to_string() );
    BGIQD::LOG::timer t(config.loger,"Overlap2Gap");

    config.LoadOverlaps();
    config.PrintCacheFreq();
    //config.FilterAllScaffItems();
    return 0;
}
