#include "soap2/contigType.h"
#include <cassert>

namespace BGIQD {
    namespace SOAP2 {

        void ContigTypeDetecter::AddContigInfo(int length, float cov )
        {
            if( length < min_len )
                return ;
            length_all += length ;
            cov_all += cov * length ;
        }

        void ContigTypeDetecter::GlobalAnalysis()
        {
            Ecov = cov_all / length_all;
            ErrorCovHigh = 1.0f ;
            HalfCovHigh = Ecov * 0.5f ;
            UniqueCovHigh = Ecov * 1.5f ;
        }

        ContigTypeDetecter::Type  ContigTypeDetecter::ContigType(int length, float cov )
        {
            if( Ecov <= 0.0f )
            {
                assert(0);
                return Type::Unknow;
            }
            if( cov < ErrorCovHigh )
            {
                return Type::Error ;
            }
            if( cov < HalfCovHigh)
            {
                return Type::Half;
            }
            if( cov < UniqueCovHigh)
            {
                return Type::Unique;
            }
            return Type::Repeat ;
        }

    }//namespace SOAP2
}//namespace BGIQD
