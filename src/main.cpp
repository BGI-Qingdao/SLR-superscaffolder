#include "contig_barcode.h"

using namespace BGIQD::JOB01;

int main()
{
    refBarcodeInfo ebi;
    refContigInfo rci;
    contigBarcodeInfo cbi;
    loadRefBarcodeInfo("../data/barcode_onref_chr19",ebi);
    loadRefContigInfo("../data/chr19.sam",rci);
    generateConrigBarcodeInfo(ebi,rci,cbi);
    printContigBarcodeInfo(cbi,"barcode_on_contig.txt");
}
