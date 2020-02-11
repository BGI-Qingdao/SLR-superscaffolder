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
#include "common/stl/mapHelper.h"
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
    std::map<unsigned int , BGIQD::stLFR::ContigBarcodeInfo> boc ;
    void LoadBarcodeOnContig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.BarcodeOnContig());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnContig to read !");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            boc[tmp.contig_id] = tmp ;
        }
        delete in;
    }
    void Init(const std::string & p )
    {
        fName.Init(p);
        barcodeOnBin.Init(1024);
        lger<<BGIQD::LOG::lstart()<<"Init finsish ..."<<BGIQD::LOG::lend();
    }


    void PrintBinSizes()
    {
        std::cout<<"ref,barcode_type,barcode_num\n";
        for(const auto & pair : boc)
        {
            auto & b2b = pair.second;
            std::set<int> count ;
            int total = 0 ;
            for( const auto & pair2 : b2b.barcodesOnPos) {
                for( int barcode : pair2.second ){
                    count.insert(barcode);
                    total++ ;
                }
            }
            std::cout<<b2b.contig_id<<','<<count.size()<<','<<total<<'\n';
        }
    }

} config;

int main(int argc ,char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix, "prefix. Input xxx.barcodeOnBin ; Output xxx.bin_cluster && xxx.cluster");
    END_PARSE_ARGS

        config.Init(prefix.to_string() );
    config.LoadBarcodeOnContig();
    config.PrintBinSizes();
    return 0;
}
