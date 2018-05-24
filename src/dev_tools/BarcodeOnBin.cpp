#include "contig_barcode.h"
#include "argsparser.h"

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main(int argc ,char **argv)
{
    initLog("BarcodeOnBin");

    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , input ,"barcodeOnContig file");
    DEFINE_ARG_REQUIRED(std::string , output, "output file");
    DEFINE_ARG_REQUIRED(int , binsize, "bin size");
    END_PARSE_ARGS

    contigBarcodeInfo cbi;
    loadContigBarcodeInfo(input.to_string() ,binsize.to_int() , cbi);
    binBarcodeInfo bbi;
    generateBinBarcodeInfo(cbi,binsize.to_int(),bbi);
    saveBinBarcodeInfo(output.to_string() , bbi);
}
