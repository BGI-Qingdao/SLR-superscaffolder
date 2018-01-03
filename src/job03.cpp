#include "contig_barcode.h"
#include "argsparser.h"

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main(int argc ,char **argv)
{
    initLog("JOB03");

    START_PARSE_ARGS
    DEFINE_ARG(std::string , input , 'i');
    DEFINE_ARG(std::string , output, 'o');
    DEFINE_ARG(int , binsize, 'b');
    END_PARSE_ARGS

    contigBarcodeInfo cbi;
    loadContigBarcodeInfo(input.to_string() ,binsize.to_int() , cbi);
    binBarcodeInfo bbi;
    generateBinBarcodeInfo(cbi,binsize.to_int(),bbi);
    saveBinBarcodeInfo(output.to_string() , bbi);
}
