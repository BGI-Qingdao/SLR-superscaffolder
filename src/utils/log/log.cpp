#include "utils/log/log.h"
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

namespace BGIQD{
namespace LOG{

    std::string logger::logstring(const std::string & str)
    {
        return module+"\t"
              +"INFO"+"\t"
              +timepoint::now().to_string() 
              +"\t:\t"+str;
    }
    timeperoid timepoint::operator-( const timepoint & prev_point ) const 
    {
        timeperoid ret ;
        ret.wall.tv_sec = wall.tv_sec - prev_point.wall.tv_sec ;
        ret.cpu = cpu - prev_point.cpu ;
        return ret ;
    }

    std::string timepoint::to_string() const 
    {
        tm * localt = localtime(&wall.tv_sec);
        //UTC   2017/12/19  15:33:33
        return std::string(tzname[0]) + "\t"
                + std::to_string(localt->tm_year+1990)+"/"
                + std::to_string(localt->tm_mon +1)+"/"
                + std::to_string(localt->tm_mday)+"\t"
                + std::to_string(localt->tm_hour)+":"
                + std::to_string(localt->tm_min )+":"
                + std::to_string(localt->tm_sec );
    }

    std::string timeperoid::to_string() const 
    {
        return std::string("wall clock : ") + std::to_string(wall.tv_sec)
               + " seconds, cpu time : "+std::to_string( cpu * 1.0f / CLOCKS_PER_SEC ) + " seconds";
    }
} //LOG
} //BGIQD
