#include "contig_barcode.h"
#include "argsparser.h"

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main(int argc ,char **argv)
{
    initLog("BarcodeOnBin");

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , input , 'i',false,"barcodeOnContig");
    DEFINE_ARG_DETAIL(std::string , output, 'o',false,"output file");
    DEFINE_ARG_DETAIL(int , binsize, 'b',false, "bin size");
    END_PARSE_ARGS

    contigBarcodeInfo cbi;
    loadContigBarcodeInfo(input.to_string() ,binsize.to_int() , cbi);
    binBarcodeInfo bbi;
    generateBinBarcodeInfo(cbi,binsize.to_int(),bbi);
    saveBinBarcodeInfo(output.to_string() , bbi);
}
