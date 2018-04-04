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

            //DEFINE_SUFFIX(

            /*
            std::string Arc() const { return m_prefix+".Arc" ;}
            std::string Arc(int round) const 
            {
                if(round == 0 )
                    return Arc() ;
                else 
                    return m_prefix+".Arc"+"_round_"+BGIQD::STRING::itos(round) ;
            }
            std::string updatedEdge() const { return m_prefix + ".updated.edge"; }
            std::string updatedEdge(int round) const 
            {
                if(round == 0 )
                    return updatedEdge() ;
                else 
                    return m_prefix+".updated.edge"+"_round_"+BGIQD::STRING::itos(round) ;
            }
            std::string contig()const { return m_prefix + ".contig" ; }
            std::string contig(int round) const 
            {
                if(round == 0 )
                    return contig() ;
                else 
                    return m_prefix+".contig"+"_round_"+BGIQD::STRING::itos(round) ;
            }
            std::string ContigIndex() const { return m_prefix + ".ContigIndex" ; }
            std::string ContigIndex(int round) const 
            {
                if(round == 0 )
                    return ContigIndex() ;
                else 
                    return m_prefix+".ContigIndex"+"_round_"+BGIQD::STRING::itos(round) ;
            }

            std::string contigroadfill() const { return m_prefix + ".contigroadfill" ; }
            std::string ContigIndex(int round) const 
            {
                if(round == 0 )
                    return ContigIndex() ;
                else 
                    return m_prefix+".ContigIndex"+"_round_"+BGIQD::STRING::itos(round) ;
            }
            std::string contigroad() const { return m_prefix + ".contigroad" ; }
            */
            private:
                std::string m_prefix;
        };
    }//namespace SOAP2
}//namespace BGIQD

#endif //__SOAP2_FILENAME_H__
