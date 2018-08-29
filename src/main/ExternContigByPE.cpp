#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "algorithm/multi_key_hash/MultiKeyHash.h"

#include "stLFR/EasySam.h"
#include "stLFR/CBB.h"

#include "soap2/graphEA.h"
#include "soap2/fileName.h"

struct AppConfig {

    BGIQD::LOG::logger loger;

    int K;

    typedef BGIQD::MultiKeyMap::BiKeyHash<unsigned int , BGIQD::EASY_SAM::PEInfo> PECache;
    PECache pe_cache ;

    BGIQD::SOAP2::GraphEA graph_ea ;
    BGIQD::SOAP2::FileNames fNames;

    std::map<unsigned int ,BGIQD::stLFR::ContigIndex> seeds;

    //PE
    void LoadPECache()
    {
        BGIQD::LOG::timer t(loger,"LoadPECache");
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.pe_pairs()) ;
        if( in == NULL )
            FATAL( "open .pe_pairs file to read failed !!! " );

        auto eachline = [this](const std::string & line) ->void 
        {
            
        };

        delete in ;
    }
    void LoadGraphEA()
    {
        BGIQD::LOG::timer t(loger,"LoadGraphEA");
        graph_ea.LoadEdge(fNames.updatedEdge(),K);
        graph_ea.LoadArc(fNames.Arc());
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.pe_seeds()) ;
        if( in == NULL )
            FATAL( "open .pe_seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            BGIQD::stLFR::ContigIndex tmp;
            tmp.InitFromString(line);
            seeds[tmp.contig] = tmp;
        }
        delete in ;
    }

    void ExternSeeds(unsigned int seeds)
    {

    }

    void ExternAll()
    {

    }

    void PrintResult()
    {

    }

    void Init( const std::string & prefix )
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
        BGIQD::stLFR::ContigIndex::K = K ;
    }
} config ;



int main(int argc , char **argv )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.pe_pairs . Output xxx.pe_contig ");
    DEFINE_ARG_REQUIRED(int , kvalue , "kvalue used by SOAP");
    END_PARSE_ARGS

    config.K = kvalue.to_int();
    config.Init(prefix.to_string());

    config.LoadSeeds();
    config.LoadGraphEA();
    config.LoadPECache();
    config.ExternAll();
    config.PrintResult();
}
