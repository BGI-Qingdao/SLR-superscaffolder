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

    refBarcodeUniqueInfo ebi;
    initLog();

    START_PARSE_ARGS
    DEFINE_ARG(std::string , input , 'i');
    DEFINE_ARG(std::string , output, 'o');
    END_PARSE_ARGS

    loadRefBarcodeUniqueInfo(input.to_string(),ebi);
    saveRefBarcodeUniqueInfo(output.to_string(),ebi);
}
