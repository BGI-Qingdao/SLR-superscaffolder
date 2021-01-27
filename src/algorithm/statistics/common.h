#ifndef __ALORITHM_STATISTICS_COMMON_H__
#define __ALORITHM_STATISTICS_COMMON_H__

#include <vector>
#include <math.h>
#include <numeric>
#include <cassert>

namespace BGIQD {
    namespace Statistics {

        template<class T, class Collection = std::vector<T> >
            void Average( const Collection & data, T & ret )
            {
                assert( data.size() > 0 );
                T total = 0 ;
                ret = std::accumulate(data.begin(),data.end(),total) / T(data.size());
            }

        template <class T, class Collection = std::vector<T> >
            void SD(const Collection & data, const T average , T & ret)
            {
                assert( data.size() > 0 );
                T sd = 0 ;
                for( const auto & d : data )
                {
                    sd += (d-average)*(d-average) ;
                }
                ret = std::sqrt( sd / T(data.size()));
            }
    }
}

#endif //__ALORITHM_STATISTICS_COMMON_H__
