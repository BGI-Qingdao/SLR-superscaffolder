#ifndef __COMMON_STL_SETHELPPER_H__
#define __COMMON_STL_SETHELPPER_H__
#include <set>
namespace BGIQD{
    namespace STL{
        template<class T>
            std::set<T> set_common( const std::set<T> & s1  , const std::set<T> &s2) {
                std::set<T> ret ;
                for( const auto & i : s1 ) 
                    if( s2.find(i) != s2.end() )
                        ret.insert(i);
                return ret ;
            }
        template<class T>
            std::set<T> set_diff_in_s1( const std::set<T> & s1  , const std::set<T> &s2) {
                std::set<T> ret ;
                for( const auto & i : s1 ) 
                    if( s2.find(i) == s2.end() )
                        ret.insert(i);
                return ret ;
            }
    }
}
#endif
