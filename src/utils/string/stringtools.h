#ifndef __COMMON_STRING_STRINGTOOLS_H__
#define __COMMON_STRING_STRINGTOOLS_H__

#include <string>
#include <vector>

/**********************************************************
 *
 * @Brief   :
 *      some useful string tools
 *          + trim
 *          + split
 *          + replace
 *
 * *******************************************************/

namespace BGIQD{
namespace STRING{

    // split by blank
    std::vector<std::string>  split(const std::string & str) ;
    std::vector<std::string>  split(const std::string & str , const std::string & spliter ) ;
    std::vector<std::string>  split(const std::string & str , const char & spliter ) ;

    // return true is all letters are digital.
    bool IsNum(const std::string & str);
} //STRING
} //BGIQD



#endif //__COMMON_STRING_STRINGTOOLS_H__
