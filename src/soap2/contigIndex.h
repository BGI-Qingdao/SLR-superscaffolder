#ifndef __SOAP2_CONTIGINDEX_H__
#define __SOAP2_CONTIGINDEX_H__

#include <string>
#include <map>
#include <iostream>

namespace BGIQD {
    namespace SOAP2 {

        struct ContigIndex
        {
            static int K ;
            unsigned int contig;
            int length ;
            bool reverse_add;

            void InitFromString(const std::string & line);
            std::string ToString() const;
        };


        struct ContigIndexMap
        {
            std::map<unsigned int , ContigIndex>  data;
            void LoadContigIndexs(std::istream & ist);
            void BuildReverseCompleteContigs();
            const ContigIndex & GetContigIndex(unsigned int id) const ;
            unsigned int BaseId(unsigned int id ) const 
            {
                return GetContigIndex(id).contig ;
            }
        };
    }
}
#endif
