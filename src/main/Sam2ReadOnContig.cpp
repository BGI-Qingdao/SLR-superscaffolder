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

    void ParseSam2ReadOnContig()
    {
        auto sam_in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig_sam());
        auto b2r_out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.read2contig());
        // basic function
        long long readId = 1 ;
        auto print1read= [&](const BGIQD::SAM::MatchData &d)
        {
            //readId\tcontigId\treadInContigPos\torigin\tbarcodeId\n
            (*b2r_out)<<readId++<<'\t'
                <<d.ref_name<<'\t'
                <<d.CalcRead1Position()<<'\t'
                <<(d.IsReverseComplete() ? '-':'+')<<'\t';
            if( has_barcode_in_read_name )
                (*b2r_out)<<BGIQD::stLFR::BarcodeIdHelper::Id( 
                        BGIQD::stLFR::readName2Barcode(d.read_name))<<'\t';
            (*b2r_out)<<(d.IsP() ? 'P' : 'E')<<'\t'
                <<(d.IsPEBothProperlyMatch() ? 'Y' : 'N')<<'\t'
                <<d.insert_size;
            (*b2r_out)<<"\n";
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
    DEFINE_ARG_OPTIONAL(std::string,barcodeList, "barcodeList file file","");
    //DEFINE_ARG_OPTIONAL(long,file_cache, "cache size of file buffer","1000000");
    END_PARSE_ARGS

    config.Init(prefix.to_string() , barcodeList.to_string() ,false );// no_stLFR.to_bool());

    BGIQD::LOG::timer t(config.loger,"Same2ReadOnContig");
    config.TryLoadBarcode2Num() ;
    config.ParseSam2ReadOnContig();
    config.PrintBarcodeList();
    return 0;
}
