#ifndef __SOAP2_CONTIGFASTA_H__
#define __SOAP2_CONTIGFASTA_H__

#include <string>
#include <map>
#include <vector>

namespace BGIQD {
    namespace SOAP2 {

        struct ContigFastA
        {
            unsigned int    id ;
            int             length ;
            float           cov ;
            int             tip;
            std::string     K;
            std::string     linear;
            int flag ;

            void Init(const std::string & line , int k) ;

            ContigFastA ReverseCompelete() const ;

            void MarkMerge()        {           flag |= 0x1 ; }
            bool IsMerge() const    { return    flag &  0x1 ; }
            void MarkSetK()         {           flag |= 0x2 ; }
            bool IsKSet() const     { return    flag &  0x2 ; }
            void MarkBase()         {           flag |= 0x4 ; }
            bool IsBase() const     { return    flag &  0x4; }
            void MarkParlindorme()  {           flag |= 0x8; }
            bool IsParlindorme()const{return    flag &  0x8; }

            void AddSeq ( const std::string & line, int k) ;

            std::string ToString() const ;
            std::string ToString(int new_id, bool marker ) const ;

            bool IsSeqComplete(int k) const 
            {
                return ( (int)K.size() == k && (int)linear.size() ==length );
            }
        }; //ConfigFasta

        class ContigFastAMap
        {
            public:

                int K;

                std::map<unsigned int , ContigFastA> contigs;

                unsigned int maxContig;

                unsigned int nextContigNum()
                {
                    maxContig += 2 ;
                    return maxContig - 1;
                }

                ContigFastA MergeContig(const std::vector<std::string> & line) ;

                ContigFastA MergeContig(const std::vector<unsigned  int> & line) ;

                void LoadContig(const std::string & file) ;


                void Init(int k )
                {
                    maxContig = 0;
                    K = k;
                }

                void buildCompeleReverse() ;
        };

    }//namespace SOAP2
}//namespace BGIQD
#endif //__SOAP2_CONTIGFASTA_H__
