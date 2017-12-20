#ifndef __LFR_UTILS_CHECK_H__
#define __LFR_UTILS_CHECK_H__

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
