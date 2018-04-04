#ifndef __SOAP2_CONTIGTYPE_H__
#define __SOAP2_CONTIGTYPE_H__

namespace BGIQD {
    namespace SOAP2 {

        struct ContigTypeDetecter
        {

            public:
            void Init(int k)
            {
                 Ecov = 0 ;
                 cov_all = 0 ;
                 length_all = 0 ;
                 K = k;
                 min_len = 2*k+1;
            }

            public:

                void AddContigInfo(int length , float cov) ;

                void GlobalAnalysis() ;

            public:

                /**
                 * Call below interface AFTER called GlobalAnalysis!!!
                 */

                enum Type
                {
                    Unknow = 0,
                    Error = 1 ,
                    Half = 2,
                    Unique = 3,
                    Repeat = 4
                };

                Type ContigType(int length, float cov );

            private:
                float Ecov;
                float cov_all;
                long length_all ;
                int K ;
                int min_len;
                float HalfCovHigh;
                float UniqueCovHigh;
                float ErrorCovHigh ;
        };

    }//namespace SOAP2
}//namespace BGIQD

#endif //__SOAP2_CONTIGTYPE_H__
