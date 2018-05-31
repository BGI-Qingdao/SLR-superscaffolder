#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/multithread/MultiThread.h"

#include "biocommon/pair/pair_sam_parser.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/barcodeId.h"
#include "stLFR/readName2Barcode.h"

#include <iostream>
#include <string>

struct AppConfig
{
    BGIQD::LOG::logger loger;
    std::string barcode_2_num_file;
    BGIQD::SOAP2::FileNames fName;
    bool has_barcode_in_read_name ;

    void Init(const std::string & prefix  ,const std::string & b2n_f, bool b )
    {
        // init loger
        BGIQD::LOG::logfilter::singleton().get("Sam2ReadInContig",BGIQD::LOG::DEBUG,loger);
        barcode_2_num_file = b2n_f ;
        fName.Init(prefix);
        has_barcode_in_read_name = ! b;
    }

    void TryLoadBarcode2Num()
    {
        if( !barcode_2_num_file.empty())
        {
            BGIQD::stLFR::BarcodeIdHelper::preload = true ;
            BGIQD::stLFR::BarcodeIdHelper::Load(barcode_2_num_file);
            loger<<BGIQD::LOG::lstart()<<" load barcodeList from "<<barcode_2_num_file<<BGIQD::LOG::lend();
        }
        else
        {
            BGIQD::stLFR::BarcodeIdHelper::preload = false ;
            loger<<BGIQD::LOG::lstart()<<" no barcodeList. will assign new barcode number ."<<BGIQD::LOG::lend();
        }
    }

    void ParseSam2ReadOnContig(std::istream & sam_in , std::ostream & b2r_out)
    {
        // basic function
        long long readId = 1 ;
        auto print1read= [&](const BGIQD::SAM::MatchData &d )
        {
            //readId\tcontigId\treadInContigPos\torigin\tbarcodeId\n
            b2r_out<<readId++<<'\t'
                <<d.ref_name<<'\t'
                <<d.first_match_position<<'\t'
                <<(d.IsReverseComplete() ? '-':'+')<<'\t';
            if( has_barcode_in_read_name )
                b2r_out<<BGIQD::stLFR::BarcodeIdHelper::Id( 
                        BGIQD::stLFR::readName2Barcode(d.read_name));
            b2r_out<<std::endl;
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
            count ++ ;
            print1read(mdata);
            if( count % 10000 == 0 )
                loger<<BGIQD::LOG::lstart()<<count<<"   pair maped reads processed ..."<<BGIQD::LOG::lend();
        };
        BGIQD::FILES::FileReaderFactory::EachLine(sam_in , parseline);
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
    DEFINE_ARG_REQUIRED(std::string,prefix, "output prefix");

    DEFINE_ARG_OPTIONAL(bool ,no_stLFR , "no barcode in read name -- if the sam file was no barcode info , open this.","");
    DEFINE_ARG_OPTIONAL(std::string,barcodeList, "barcodeList file file","");
    DEFINE_ARG_OPTIONAL(long,file_cache, "cache size of file buffer","1000000");
    END_PARSE_ARGS

    if( file_cache.setted && file_cache.to_int() > 1024)
    {
        BGIQD::FILES::FileReaderFactory::ResizeBuff(std::cin , file_cache.to_int());
        BGIQD::FILES::FileWriterFactory::ResizeBuff(std::cout , file_cache.to_int());
    }

    config.Init(prefix.to_string() , barcodeList.to_string() , no_stLFR.to_bool());

    BGIQD::LOG::timer t(config.loger,"Same2ReadOnContig");
    config.TryLoadBarcode2Num() ;
    config.ParseSam2ReadOnContig(std::cin,std::cout);
    config.PrintBarcodeList();
    return 0;
}
