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

            // SOAP2
            DEFINE_SUFFIX(Arc,".Arc");

            DEFINE_SUFFIX(updatedEdge,".updated.edge");

            DEFINE_SUFFIX(contig, ".contig");

            DEFINE_SUFFIX(ContigIndex, ".ContigIndex");
            // StaticUnique
            DEFINE_SUFFIX(seeds, ".seeds");
            // Same2ReadOnContig
            DEFINE_SUFFIX(read2contig_sam, ".read2contig.sam");

            DEFINE_SUFFIX(read2contig, ".read2contig");
            // ChopBin
            DEFINE_SUFFIX(barcodeList, ".barcodeList");

            DEFINE_SUFFIX(BarcodeOnBin, ".barcodeOnBin");

            DEFINE_SUFFIX(BarcodeOnContig, ".barcodeOnContig");
            // Bin Cluster
            DEFINE_SUFFIX(cluster, ".cluster");

            DEFINE_SUFFIX(bin_cluster, ".bin_cluster");
            // ContigDLink
            DEFINE_SUFFIX(connInfo, ".connInfo");
            // LinearCDG
            DEFINE_SUFFIX(contigroad , ".contigroad");
            // FillContigRoad
            DEFINE_SUFFIX(contigroadfill,".contigroadfill");
            // MergeContig
            DEFINE_SUFFIX(super_used, ".super_used");

            DEFINE_SUFFIX(super_only, ".super_only");

            DEFINE_SUFFIX(super_and_left, ".super_and_left");
            private:
                std::string m_prefix;
        };
    }//namespace SOAP2
}//namespace BGIQD

#endif //__SOAP2_FILENAME_H__
