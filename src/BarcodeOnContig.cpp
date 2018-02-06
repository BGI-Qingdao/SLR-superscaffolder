#include "contig_barcode.h"
#include "argsparser.h"

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

int main(int argc , char ** argv)
{
    initLog("BarcodeOnContig");

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , refBarcode , 'i',false,"barcodeOnRef");
    DEFINE_ARG_DETAIL(std::string , refContig, 'c',false, "contigOnRef ");
    DEFINE_ARG_DETAIL(std::string , output, 'o',false,"output prefix");
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
