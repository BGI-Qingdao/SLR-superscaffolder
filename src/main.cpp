#include "contig_barcode.h"

using namespace BGIQD::JOB01;

int main()
{
    refBarcodeInfo ebi;
    refContigInfo rci;
    contigBarcodeInfo cbi;
    loadRefBarcodeInfo("",ebi);
    loadRefContigInfo("",rci);
    generateConrigBarcodeInfo(ebi,rci,cbi);
    printContigBarcodeInfo(cbi,"");
}
