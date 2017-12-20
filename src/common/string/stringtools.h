#ifndef __COMMON_STRING_STRINGTOOLS_H__
#define __COMMON_STRING_STRINGTOOLS_H__

#include <string>
#include <vector>

namespace BGIQD{
namespace STRING{

std::vector<std::string>  split(const std::string & str , const std::string & spliter ) ;

std::string ltrim(const std::string & str);
std::string rtrim(const std::string & str);
std::string trim(const std::string & str);
} //STRING
} //BGIQD



#endif //__COMMON_STRING_STRINGTOOLS_H__
