#ifndef __SOAP2_FILENAME_H__
#define __SOAP2_FILENAME_H__

#include <string>
#include "common/string/stringtools.h"

namespace BGIQD {
    namespace SOAP2 {

        struct FileNames
        {
            void Init(const std::string & prefix)
            {
                m_prefix = prefix ;
            }

#define DEFINE_SUFFIX(name,suffix) \
            std::string name() const { return m_prefix + suffix ; } \
            std::string name(int round) const \
            {\
                if ( round == 0 )\
                {\
                    return name();\
                }\
                else\
                {\
                    return m_prefix + suffix +"_round_"+BGIQD::STRING::itos(round) ;\
                }\
            }

            DEFINE_SUFFIX(Arc,".Arc");

            DEFINE_SUFFIX(updatedEdge,".updated.edge");

            DEFINE_SUFFIX(contig, ".contig");

            DEFINE_SUFFIX(ContigIndex, ".ContigIndex");

            DEFINE_SUFFIX(contigroadfill,".contigroadfill");

            DEFINE_SUFFIX(contigroad , ".contigroad");

            DEFINE_SUFFIX(super_used, ".super_used");

            DEFINE_SUFFIX(super_only, ".super_only");

            DEFINE_SUFFIX(super_and_left, ".super_and_left");

            DEFINE_SUFFIX(barcodeList, ".barcodeList");

            DEFINE_SUFFIX(BarcodeOnContig, ".barcodeOnContig");

            DEFINE_SUFFIX(BarcodeOnBin, ".barcodeOnBin");

            DEFINE_SUFFIX(seeds, ".seeds");

            private:
                std::string m_prefix;
        };
    }//namespace SOAP2
}//namespace BGIQD

#endif //__SOAP2_FILENAME_H__
