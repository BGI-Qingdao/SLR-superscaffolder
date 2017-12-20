#include "stringtools.h"

namespace BGIQD{
namespace STRING{

    std::vector<std::string>  split( const std::string & str , const std::string & spliter ) 
    {
        std::vector<std::string> ret;
        std::string::size_type pos1 = 0;
        std::string::size_type pos2 = str.find(spliter);
        while( pos2 != str.npos )
        {
            ret.push_back( str.substr( pos1 , pos2-pos1 ) );
            pos1 = pos2 + spliter.length() ;
            pos2 = str.find( spliter , pos1 ) ;
        }
        if( pos1 != str.length() )
        {
            ret.push_back(str.substr(pos1));
        }
        return ret ;
    }

}
}
