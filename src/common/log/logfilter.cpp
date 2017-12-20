#include "logfilter.h"
#include <iostream>

namespace BGIQD{
namespace LOG {

    logfilter* logfilter::the_one = NULL ;
    //TODO : not thread safe here.
    logfilter & logfilter::singleton() 
    {
        if( the_one == NULL )
        {
            the_one = new logfilter();
        }
        return *the_one;
    }

    logger * logfilter::get(const std::string module, BGIQD::LOG::loglevel level )
    {
        themes.push_back((ilogtheme*)( new detaillog(module,levelname(level))));
        //TODO : parse argument and filter!
        return new logger(&std::cerr, themes.back());
    }

    logfilter::~logfilter()
    {
        for(size_t i = 0 ; i < themes.size() ; i++)
        {
            delete themes[i];
            themes[i] = NULL;
        }
    }
}
}
