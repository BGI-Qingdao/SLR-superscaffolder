/**********************************************************
 *
 * @Brief  :
 *
 *   Parse SAM and print in customized format.
 *
 * *******************************************************/
#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"
#include "utils/multithread/MultiThread.h"
#include "utils/misc/contigIndex.h"
#include "utils/misc/fileName.h"
#include "utils/misc/TagId.h"
#include "utils/sam/EasySam.h"
#include "utils/sam/sam_parser.h"
#include "utils/seq/fastq.h"

#include <iostream>
#include <string>
#include <cassert>

//
// Struct to wrap all global variables and functions
//
struct AppConfig
{
    BGIQD::LOG::logger loger;
    std::string barcode_2_num_file;
    BGIQD::MISC::FileNames fName;
    bool has_barcode_in_read_name ;

    typedef BGIQD::SEQ::stLFRHeader stLFRHeader ;

    void Init(const std::string & prefix  ,const std::string & b2n_f, bool b )
    {
        // init loger
        BGIQD::LOG::logfilter::singleton().get("Sam2ReadInContig",BGIQD::LOG::DEBUG,loger);
        barcode_2_num_file = b2n_f ;
        fName.Init(prefix);
        has_barcode_in_read_name = ! b;
    }

    BGIQD::MISC::StringIdCache barcodeIds;

    BGIQD::MISC::StringIdCache readNameIds ;

    void LoadBarcode2Num()
    {
        BGIQD::LOG::timer t (loger,"LoadBarcode2Num");
        barcodeIds.preload = true ;
        barcodeIds.Load(fName.barcodeList());
    }

    void LoadRead2Num()
    {
        BGIQD::LOG::timer t(loger,"LoadRead2Num");
        readNameIds.preload = true ;
        readNameIds.Load(fName.readNameList());
    }

    void ParseSam2ReadOnContig()
    {
        std::vector<BGIQD::EASY_SAM::EasySam> easy_cache;
        auto sam_in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig_sam());
        // basic function
        // long long readId = 1 ;
        auto print1read= [&](const BGIQD::SAM::MatchData &d)
        {
            BGIQD::EASY_SAM::EasySam tmp;
            //tmp.read_id = readId++;
            tmp.contig_name = std::stoul(d.ref_name);
            tmp.pos_1bp = d.CalcRead1Position();
            tmp.match_reverse = d.IsReverseComplete() ;
            stLFRHeader header;
            header.Init(d.read_name);
            assert(header.type != stLFRHeader::ReadType::Unknow );
            tmp.barcode = barcodeIds.Id(header.barcode_str);
            tmp.read_id = readNameIds.Id(header.readName);
            tmp.is_p = d.IsP();
            if( tmp.is_p == 0 )
                tmp.read_id ++ ;
            tmp.pe_match = (d.IsPEBothMatch() && ! d.XA);
            tmp.insert_size = d.insert_size ;
            easy_cache.push_back(tmp);
        };
        long long count = 0 ;
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
            if( ! mdata.IsPrimaryMatch() )
                return ;
            count ++ ;
            print1read(mdata);
            if( count % 1000000 == 0 )
                loger<<BGIQD::LOG::lstart()<<count<<"   pair maped reads processed ..."<<BGIQD::LOG::lend();
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*sam_in , parseline);
        delete sam_in ;

        auto b2r_out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.read2contig());
        for( const auto & item : easy_cache)
        {
            (*b2r_out)<<item.ToString()<<'\n';
        }
        delete b2r_out;
    }
}config;

int main(int argc , char ** argv)
{
    // parse args
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string,prefix, "prefix. Input xxx.read2contig.sam \n\
                                                                              xxx.readNameList\n\
                                                                              xxx.barcodeList ;\n\
                                                                    Output xxx.read2contig");
    END_PARSE_ARGS

    config.Init(prefix.to_string() , ""/*barcodeList.to_string()*/ ,false );// no_stLFR.to_bool());
    BGIQD::LOG::timer t(config.loger,"Same2ReadOnContig");
    config.LoadRead2Num();
    config.LoadBarcode2Num() ;
    config.ParseSam2ReadOnContig();
    return 0;
}
