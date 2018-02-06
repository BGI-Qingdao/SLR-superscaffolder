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
    DEFINE_ARG_DETAIL(std::string , input , 'i',false,"barcodeOnRef_wenchao");
    DEFINE_ARG_DETAIL(std::string , output, 'o',false,"output");
    END_PARSE_ARGS

    refBarcodeUniqueInfo ebi;

    loadRefBarcodeUniqueInfo(input.to_string(),ebi);
    saveRefBarcodeUniqueInfo(output.to_string(),ebi);
}
