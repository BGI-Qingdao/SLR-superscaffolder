#ifndef __ALGORITHM_DISTRIBUTION_DISTRIBUTION_H__
#define __ALGORITHM_DISTRIBUTION_DISTRIBUTION_H__
#include <map>
#include "algorithm/interval/Interval.h"
namespace BGIQD {
    namespace DISTRIBUTION {

        template<class T>
            struct IntervalDistribution
            {
                typedef BGIQD::INTERVAL::Interval<T>  Item;
                std::map<Item, int> freqs;
                std::map<Item, float> percents;

                void Init(const T & bin ,const  T &  min ,const T & max )
                {
                    for( T m = min ; m < max ; m += bin )
                    {
                        freqs[Item(m,m+bin-1)] = 0 ;
                    }
                }

                void Count(const T & t , int count  =1 )
                {
                    for( auto & pair : freqs )
                    {
                        if ( pair.first.IsContain(t) )
                        {
                            pair.second += count ;
                            return ;
                        }
                    }
                }

                void CalcPercent()
                {
                    int total = 0 ;
                    for( auto & pair : freqs )
                    {
                        total += pair.second ;
                    }

                    for( auto & pair : freqs )
                    {
                        percents[pair.first] = (float)pair.second  / (float)total;
                    }
                }

                float GetPercent(const T & t)
                {
                    for( auto & pair : percents)
                    {
                        if ( pair.first.IsContain(t) )
                        {
                            return pair.second ;
                        }
                    }
                    return 0.0f ;
                }

                std::string ToString()
                {
                    std::string ret ;
                    for( auto & pair : percents)
                    {
                        ret += pair.first.ToString() ;
                        ret += '\t';
                        ret += std::to_string(pair.second);
                        ret += '\n';
                    }
                    return ret ;
                }
            };
    }

}

#endif //__ALGORITHM_DISTRIBUTION_DISTRIBUTION_H__
