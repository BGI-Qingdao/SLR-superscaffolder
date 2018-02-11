#ifndef __COMMON_FREQ_FREQ_H__
#define __COMMON_FREQ_FREQ_H__
#include <map>
#include <sstream>
namespace BGIQD{
namespace FREQ{

    template< class Key >
        class Freq
        {
            public:
                void Touch( const Key k)
                {
                    if( data.find(k) != data.end() )
                        data[k] ++ ;
                    else
                        data[k] = 1 ;
                }


                std::string ToString()
                {
                    std::ostringstream ret ;
                    for( const auto i : data )
                        ret<<i.first<<'\t'<<i.second<<std::endl;
                    return ret.str();
                }
                std::map<Key, long> data;
        };
}
}


#endif
