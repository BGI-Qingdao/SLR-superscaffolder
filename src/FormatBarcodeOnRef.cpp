#include "contig_barcode.h"
#include "argsparser.h"


using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;
/*****************************************************************************
 *
 * remove PCR duplcate from barcode_on_ref_chr19
 *
 * ***************************************************************************/

int main(int argc , char **argv )
{
    initLog("FormatBarcodeOnRef");

    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , input ,"barcodeOnRef_wenchao");
    DEFINE_ARG_REQUIRED(std::string , output,"output");
    END_PARSE_ARGS

    refBarcodeUniqueInfo ebi;

    loadRefBarcodeUniqueInfo(input.to_string(),ebi);
    saveRefBarcodeUniqueInfo(output.to_string(),ebi);
}
