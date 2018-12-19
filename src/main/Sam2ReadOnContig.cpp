#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/multithread/MultiThread.h"

#include "biocommon/pair/pair_sam_parser.h"
#include "biocommon/fastq/fastq.h"


#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/barcodeId.h"
#include "stLFR/readName2Barcode.h"
#include "stLFR/EasySam.h"
#include "stLFR/StringIdCache.h"

#include <iostream>
#include <string>
#include <cassert>

struct AppConfig
{
    BGIQD::LOG::logger loger;
    std::string barcode_2_num_file;
    BGIQD::SOAP2::FileNames fName;
    bool has_barcode_in_read_name ;

    typedef BGIQD::FASTQ::stLFRHeader stLFRHeader ;

    void Init(const std::string & prefix  ,const std::string & b2n_f, bool b )
    {
        // init loger
        BGIQD::LOG::logfilter::singleton().get("Sam2ReadInContig",BGIQD::LOG::DEBUG,loger);
        barcode_2_num_file = b2n_f ;
        fName.Init(prefix);
        has_barcode_in_read_name = ! b;
    }

    BGIQD::stLFR::StringIdCache barcodeIds;

    BGIQD::stLFR::StringIdCache readNameIds ;

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
            //tmp.barcode = BGIQD::stLFR::BarcodeIdHelper::Id(BGIQD::stLFR::readName2Barcode(d.read_name));
            tmp.barcode = barcodeIds.Id(header.barcode_str);
            tmp.read_id = readNameIds.Id(header.readName);
            tmp.is_p = d.IsP();
            if( tmp.is_p == 0 )
                tmp.read_id ++ ;
            tmp.pe_match = (d.IsPEBothMatch() && ! d.XA);
            tmp.insert_size = d.insert_size ;
            easy_cache.push_back(tmp);
        };
        // parse sam and print
        //BGIQD::SAM::PairedSAMParser parser(sam_in);
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

    void PrintBarcodeList()
    {
        // print barcodeList
        if( barcode_2_num_file.empty() )
        {
            BGIQD::stLFR::BarcodeIdHelper::Print(fName.barcodeList());
        }
    }

}config;



int main(int argc , char ** argv)
{
    // parse args
    START_PARSE_ARGS
    //DEFINE_ARG_OPTIONAL(bool ,no_stLFR , "no barcode in read name -- if the sam file was no barcode info , open this.","");
    DEFINE_ARG_REQUIRED(std::string,prefix, "prefix. Input xxx.read2contig.sam ; Output xxx.read2contig && xxx.barcodeList");
    //DEFINE_ARG_OPTIONAL(std::string,barcodeList, "barcodeList file file","");
    //DEFINE_ARG_OPTIONAL(long,file_cache, "cache size of file buffer","1000000");
    END_PARSE_ARGS

    config.Init(prefix.to_string() , ""/*barcodeList.to_string()*/ ,false );// no_stLFR.to_bool());
    BGIQD::LOG::timer t(config.loger,"Same2ReadOnContig");
    config.LoadRead2Num();
    config.LoadBarcode2Num() ;
    config.ParseSam2ReadOnContig();
    //config.PrintBarcodeList();
    return 0;
}
