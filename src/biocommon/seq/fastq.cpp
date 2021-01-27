#include "biocommon/seq/fastq.h"
#include "common/string/stringtools.h"

#include <string>
#include <sstream>

namespace BGIQD {

    namespace SEQ {

        void stLFRHeader::Init( const std::string & line ) 
        {
            Reset() ;

            auto items = BGIQD::STRING::split(line ,'@');

            assert(!items.empty());

            auto items_1 = BGIQD::STRING::split(items[0]);

            assert(!items_1.empty());

            auto items_2 = BGIQD::STRING::split(items_1[0],'#');

            assert(items_2.size() == 2);

            readName = items_2[0] ;

            auto items_3 = BGIQD::STRING::split(items_2[1] , '/');

            assert( !items_3.empty());

            barcode_str = items_3[0] ;

            type = readName_barcodeStr ;

            if( items_3.size() == 2 && BGIQD::STRING::IsNum(items_3[1]) )
            {
                readIndex = std::stoi(items_3[1]);
                type = readName_barcodeStr_index ;
            }

            if( items_1.size() >1 && BGIQD::STRING::IsNum(items_1[1]) )
            {
                barcode_num = std::stoi(items_1[1]);
                type = readName_barcodeStr_index_barcodeNum ;
            }
        }

        std::string stLFRHeader::Head() const 
        {
            assert(type != Unknow );
            std::ostringstream ost ;
            ost<<'@'<<readName<<'#'<<barcode_str;
            if( type == readName_barcodeStr_index )
                ost<<'/'<<readIndex;
            if( type == readName_barcodeStr_index_barcodeNum )
                ost<<'/'<<readIndex<<'\t'<<barcode_num;
            return ost.str();
        }

    }
}
