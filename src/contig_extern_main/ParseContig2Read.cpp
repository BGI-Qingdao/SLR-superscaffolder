#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/multithread/MultiThread.h"
#include "common/error/Error.h"

#include "biocommon/sam_bam/sam_parser.h"
#include "biocommon/fastq/fastq.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/StringIdCache.h"
#include "stLFR/readName2Barcode.h"
#include "stLFR/EasySam.h"

#include "algorithm/incr_array/incr_array.h"

#include <iostream>
#include <string>

struct AppConfig
{
    typedef std::vector<BGIQD::EASY_SAM::EasySam_V1> EasySamCache;

    typedef BGIQD::INCRARRAY::IncrArray<EasySamCache> PBuffer;

    PBuffer print_buffer;

    BGIQD::LOG::logger loger;

    BGIQD::SOAP2::FileNames fName;

    BGIQD::stLFR::StringIdCache barcodeIds;

    BGIQD::stLFR::StringIdCache readNameIds ;

    BGIQD::MultiThread::MultiThread mt;

    void Init(const std::string & prefix)
    {
        // init loger
        BGIQD::LOG::logfilter::singleton().get("Sam2ReadInContig",BGIQD::LOG::DEBUG,loger);
        fName.Init(prefix);
        print_buffer.Init(10000);
    }

    void LoadBarcode2Num()
    {
        BGIQD::LOG::timer t(loger,"LoadBarcode2Num");
        barcodeIds.preload = true ;
        barcodeIds.Load(fName.barcodeList());
    }

    void LoadRead2Num()
    {
        BGIQD::LOG::timer t(loger,"LoadRead2Num");
        readNameIds.preload = true ;
        readNameIds.Load(fName.readNameList());
    }

    void ParseContig2read_Sam(const std::string & file )
    {
        int buffer_index = 0 ;
        int cache_size = 30000 ;
        long long count = 0 ;
        EasySamCache easy_cache;

        // basic function 1 : save data to print buffer
        auto save_buffer = [&] () {
            print_buffer.push_back(EasySamCache());
            auto & buffer = print_buffer[buffer_index];
            std::swap( buffer , easy_cache);
            std::function<void()> job = std::bind(
                    &AppConfig::PrintEasySam 
                    , this
                    , buffer_index);
            buffer_index ++ ;
            mt.AddJob(job);
        };

        // basic function 2 : make a MatchData to EasySam
        auto print1read= [&](const BGIQD::SAM::MatchData &d)
        {
            BGIQD::EASY_SAM::EasySam_V1 tmp;

            BGIQD::FASTQ::stLFRHeader header;

            //header.Init(d.read_name);
            header.Init(d.ref_name);

            tmp.read_id = readNameIds.Id(header.readName);
            if( header.type == 
                    BGIQD::FASTQ::stLFRHeader::
                    ReadType::readName_barcodeStr_index
                    || header.type == 
                    BGIQD::FASTQ::stLFRHeader::
                    ReadType::readName_barcodeStr_index_barcodeNum
              )
            {
                tmp.read_index = header.readIndex;
            }
            else
            {
                tmp.read_index = d.IsP() ? 1 : 2 ;
            }
            tmp.contig_name = std::stoul(d.read_name);
            // TODO : this is the contig's left most position in read !!!
            tmp.left_1bp= d.CalcLeft1Position();
            tmp.match_reverse = d.IsReverseComplete() ;
            tmp.barcode = barcodeIds.Id(header.barcode_str);

            easy_cache.push_back(tmp);

            if( int(easy_cache.size()) == cache_size)
            {
                save_buffer();
            }
        };
        // basic function 3 : parse sam to MatchData
        auto parseline = [&](const std::string & line)
        {
            BGIQD::SAM::LineParser l(line);
            if( ! l.IsVaid() || l.IsHead() )
            {
                return ;
            }
            auto mdata = l.ParseAsMatchData();
            if( mdata.UnMap())
                return ;
            if(  mdata.IsSupplementaryMatch() )
                return ;
            count ++ ;
            print1read(mdata);
            if( count % 1000000 == 0 )
                loger<<BGIQD::LOG::lstart()<<count<<"   pair maped reads processed ..."<<BGIQD::LOG::lend();
        };
        // function logic :
        auto sam_in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(file);
        if( sam_in == NULL )
            FATAL(" failed to open xxx.contig2rx.sam to read !!!");
        BGIQD::FILES::FileReaderFactory::EachLine(*sam_in , parseline);
        delete sam_in ;
        if( easy_cache.size() > 0 )
        {
            save_buffer();
        }
    }

    void ParseContig2r1_Sam()
    {
        BGIQD::LOG::timer t(loger,"Contig2r1");
        ParseContig2read_Sam(fName.contig2r1_sam());
    }

    void ParseContig2r2_Sam()
    {
        BGIQD::LOG::timer t(loger,"Contig2r2");
        ParseContig2read_Sam(fName.contig2r2_sam());
    }

    std::ostream * out ;
    void StartWriteThread()
    {
        out = NULL ;
        out = BGIQD::FILES::
              FileWriterFactory::GenerateWriterFromFileName(
                      fName.contig2read_v1());
        if( out == NULL )
            FATAL("failed to open xxx.read2contig_v1 to write");
        mt.Start(1);
    }

    void EndWriteThread()
    {
        mt.End();
        mt.WaitingStop();
        delete out ;
        out = NULL ;
    }

    void PrintEasySam( int i )
    {
        auto & buffer = print_buffer[i] ;
        for( const auto & i : buffer)
        {
            (*out)<<i.ToString()<<'\n';
        }
        buffer.clear() ;
        buffer.shrink_to_fit();
    }

}config;

int main(int argc , char ** argv)
{
    // parse args
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,prefix, "prefix. \n\
                                                Input\n\
                                                    xxx.contig2r1.sam\n\
                                                    xxx.contig2r2.sam\n\
                                                    xxx.barcodeList\n\
                                                    xxx.readNameList\n\
                                                Output xxx.contig2read_v1");
    END_PARSE_ARGS

    config.Init(prefix.to_string());

    BGIQD::LOG::timer t(config.loger,"ParseContig2eRead");

    config.LoadBarcode2Num() ;
    config.LoadRead2Num();

    config.StartWriteThread();
    config.ParseContig2r1_Sam();
    config.ParseContig2r2_Sam();
    config.EndWriteThread();
    return 0;
}
