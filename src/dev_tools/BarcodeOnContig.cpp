#include "contig_barcode.h"
#include "common/args/argsparser.h"

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

int main(int argc , char ** argv)
{
    initLog("BarcodeOnContig");

    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , refBarcode , "barcodeOnRef");
    DEFINE_ARG_REQUIRED(std::string , refContig, "contigOnRef ");
    DEFINE_ARG_REQUIRED(std::string , output, "output prefix");
    END_PARSE_ARGS

    refBarcodeInfo ebi;
    refContigInfo rci;
    contigBarcodeInfo cbi;
    contigLens lens;
    loadRefBarcodeInfo(refBarcode.to_string(),ebi);
    loadRefContigInfo(refContig.to_string(),rci,lens);
    generateConrigBarcodeInfo(ebi,rci,cbi);
    printContigBarcodeInfo(cbi,lens,output.to_string());
}
