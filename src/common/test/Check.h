#ifndef __LFR_UTILS_CHECK_H__
#define __LFR_UTILS_CHECK_H__
#include "log.h"


#define CHECK_STRUCT( expect , value )\
    if( expect != value ) \
        std::cerr<<__FILE__<<":"<<__LINE__<<" not match !!"<<std::endl; 

#define CHECK_STRUCT_AND_ONERR( expect , value , action) \
    if( expect != value ) \
    {\
        std::cerr<<__FILE__<<":"<<__LINE__<<" not match !!"<<std::endl;\
        action \
    }

#define CHECK( expect , value ) \
    if( expect != value ) \
        std::cerr<<__FILE__<<":"<<__LINE__<<" expect "<<expect<<" but "<<value<<std::endl; 

#define CHECK_AND_ONERR( expect , value , action) \
    if( expect != value ) \
    {\
        std::cerr<<__FILE__<<":"<<__LINE__<<" expect "<<expect<<" but "<<value<<std::endl;\
        action \
    }


#endif //__LFR_UTILS_CHECK_H__
