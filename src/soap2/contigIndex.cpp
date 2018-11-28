#include "soap2/contigIndex.h"
#include <sstream>
namespace BGIQD {
    namespace SOAP2 {

        int ContigIndex::K  ;

        std::string ContigIndex::ToString() const 
        {
            std::ostringstream ist;
            ist<<contig<<'\t'<<length<<'\t'<<reverse_add;
            return ist.str();
        }
        void ContigIndex::InitFromString(const std::string &line )
        {
            std::istringstream ist(line);
            ist>>contig>>length;
            if( ! ist.eof() )
                ist>>reverse_add;
            else
                reverse_add = 1 ;
        }

    }
}
