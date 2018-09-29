#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"
#include "common/error/Error.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "soap2/fileName.h"
#include "stLFR/CBB.h"
#include "stLFR/EasySam.h"

#include <vector>
#include <map>

struct AppConfig
{
    BGIQD::SOAP2::FileNames fName ;
    BGIQD::LOG::logger log;

    void Init( const std::string & prefix )
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("FakePE2BC",BGIQD::LOG::loglevel::INFO,log);
    }

    std::map<unsigned int ,BGIQD::stLFR::ContigIndex> seeds;
    std::map<unsigned int , BGIQD::stLFR::ContigBarcodeInfo> boc;

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(log,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            BGIQD::stLFR::ContigIndex tmp ;
            tmp.InitFromString(line) ;
            seeds[tmp.contig] = tmp ;
            boc[tmp.contig].contig_id= tmp.contig;
        }
        delete in ;
    }

    
    void LoadPEPair()
    {
        BGIQD::LOG::timer t(log,"LoadPECahce");
        BGIQD::FREQ::Freq<int> ISFreq;
        std::string line ;

        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.pe_pairs());
        if ( in == NULL )
            FATAL(" open xxx.pe_pair to read failed !!! ");

        while( ! std::getline( *in , line ).eof() )
        {
            BGIQD::EASY_SAM::PEInfo tmp ;
            tmp.InitFromString(line);
            if( boc.find(tmp.contig1) == boc.end()   )
                continue ;
            if( boc.find(tmp.contig2) == boc.end()   )
                continue ;
            // Use read1 represent this read pair
            boc[tmp.contig1].Touch(tmp.pos_1bp1,tmp.read1);
            boc[tmp.contig2].Touch(tmp.pos_1bp2,tmp.read1);
        }
    }
    void PrintFakeBarcodeOnContig()
    {
        auto in1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.BarcodeOnContig_fake());
        for( const auto & pair : boc)
        {
            (*in1)<<pair.second.ToString()<<'\n';
        }
        delete in1 ;
    }

} config ;



int  main(int argc, char **argv)
{
    //step0 Parse parmeters...
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix of filenames .\n\
                                                    In \n\
                                                        xxx.pe_pair;\n\
                                                        xxx.seeds\n\
                                                    Output \n\
                                                        xxx.barcodeOnContig_fake");
    END_PARSE_ARGS;

    config.Init( prefix.to_string());
    config.LoadSeeds() ;
    config.LoadPEPair();
    config.PrintFakeBarcodeOnContig();
}
