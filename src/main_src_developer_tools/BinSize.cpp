#include "algorithm/incr_array/incr_array.h"
#include "algorithm/collection/collection.h"

#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/multithread/MultiThread.h"
#include "common/stl/mapHelper.h"
#include "common/log/log.h"
#include "common/error/Error.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/middle_valiad/MiddleValid.h"
#include "common/freq/freq.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <algorithm>
#include <iostream>
#include <string>

struct AppConfig
{

    BGIQD::stLFR::BarcodeOnBinArray barcodeOnBin ;
    BGIQD::SOAP2::FileNames fName;
    BGIQD::LOG::logger lger;

    void Init(const std::string & p )
    {
        fName.Init(p);
        barcodeOnBin.Init(1024);
        lger<<BGIQD::LOG::lstart()<<"Init finsish ..."<<BGIQD::LOG::lend();
    }

    void LoadB2BArray( )
    {
        BGIQD::stLFR::LoadBarcodeOnBinArray( fName.BarcodeOnBin(middle_name) , barcodeOnBin);
    }

    BGIQD::FREQ::Freq<int> type_Freq;
    void PrintBinSizes()
    {

        std::cout<<"ref,bin_id,barcode_type\n";
        for(size_t i = 0 ; i <barcodeOnBin.size() ; i++)
        {
            auto & b2b = barcodeOnBin.at(i);
            std::cout<<b2b.contigId<<','<<b2b.binId<<','<<b2b.collections.keysize()<<'\n';
        }
    }
    std::string middle_name ;

} config;

int main(int argc ,char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix. Input xxx.barcodeOnBin ; Output xxx.bin_cluster && xxx.cluster");
    DEFINE_ARG_OPTIONAL(std::string,  middle_name, "the middle name of output suffix " ,"");
    END_PARSE_ARGS

    config.middle_name = middle_name.to_string() ;
    config.Init(prefix.to_string() );
    config.LoadB2BArray() ;
    config.PrintBinSizes();
    return 0;
}
