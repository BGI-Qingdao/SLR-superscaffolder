#include "contig_barcode.h"
#include "argsparser.h"

using namespace BGIQD::ARGS;
using namespace BGIQD::JOB01;

int main(int argc , char ** argv)
{
    initLog("JOB01");

    START_PARSE_ARGS
    DEFINE_ARG(std::string , refBarcode , 'i');
    DEFINE_ARG(std::string , refContig, 'c');
    DEFINE_ARG(std::string , output, 'o');
    END_PARSE_ARGS

    refBarcodeInfo ebi;
    refContigInfo rci;
    contigBarcodeInfo cbi;
    loadRefBarcodeInfo(refBarcode.to_string(),ebi);
    loadRefContigInfo(refContig.to_string(),rci);
    generateConrigBarcodeInfo(ebi,rci,cbi);
    printContigBarcodeInfo(cbi,output.to_string());
}
