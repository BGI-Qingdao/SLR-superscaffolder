#include "biocommon/pair/pair_sam_parser.h"
#include "stLFR/barcodeId.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "stLFR/readName2Barcode.h"
#include <iostream>

static BGIQD::LOG::logger loger;
//using namespace BGIQD::ARGS;
int main(int argc , char ** argv)
{

    START_PARSE_ARGS
    DEFINE_ARG(bool ,no_stLFR , 'n');
    DEFINE_ARG(std::string,barcodeList, 'b');
    DEFINE_ARG(std::string,prefix, 'o');
    END_PARSE_ARGS

    if( barcodeList.setted )
    {
        BGIQD::stLFR::BarcodeIdHelper::preload = true ;
        BGIQD::stLFR::BarcodeIdHelper::Load(barcodeList.to_string());
        loger<<BGIQD::LOG::lstart()<<" load barcodeList from "<<barcodeList.to_string()<<BGIQD::LOG::lend();
    }
    else
    {
        BGIQD::stLFR::BarcodeIdHelper::preload = false ;
        loger<<BGIQD::LOG::lstart()<<" no barcodeList. will assign new barcode number ."<<BGIQD::LOG::lend();
    }

    BGIQD::LOG::logfilter::singleton().get("Sam2ReadInContig",BGIQD::LOG::DEBUG,loger);
    BGIQD::SAM::PairedSAMParser parser(std::cin);
    long long readId = 1 ;

    auto print1read= [&](const BGIQD::SAM::MatchData &d )
    {
        //readId\tcontigId\treadInContigPos\torigin\tbarcodeId\n
        std::cout<<readId++<<'\t'
            <<d.ref_name<<'\t'
            <<d.first_match_position<<'\t'
            <<(d.IsReverseComplete() ? '-':'+')<<'\t';
        if( ! no_stLFR.to_bool() )
            std::cout<<BGIQD::stLFR::BarcodeIdHelper::Id( 
                    BGIQD::stLFR::readName2Barcode(d.read_name));
        std::cout<<std::endl;
    };
    long long count = 0 ;
    while(1)
    {
        auto p = parser.CurrentPair();
        if( ! p.first.Valid() || !p.second.Valid() )
            break;
        if( p.first.UnMap() || p.second.UnMap() )
            continue;
        count ++ ;
        print1read(p.first);
        print1read(p.second);
        if( count % 1000 == 0 )
            loger<<BGIQD::LOG::lstart()<<count<<"   pair maped reads processed ..."<<BGIQD::LOG::lend();
    }
    if( ! barcodeList.setted && prefix.setted )
    {
        BGIQD::stLFR::BarcodeIdHelper::Print(prefix.to_string() + ".barcodeList");
    }
    return 0;
}
