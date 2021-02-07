#ifndef __LFR_TEST_TEST_H__
#define __LFR_TEST_TEST_H__

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <functional>

#include "utils/log/log.h"

/**********************************************************
 *
 * A simple framework for unittest
 *
 * How to use it ? 
 *
 *   + the main file is in main/test.cpp
 *   + write a xxx_test.cpp and link it after compiler then all test code will automatic hooked into bin/test
 *   + find examples in other folders under utils to add more unittest.
 *
 * A simple example
 *
 *   TEST_MODULE(example)
 *   TEST(test_01)
 *   {
 *       CHECK(1 ,3-2)   <-- this check will pass 
 *       CHECK(2 ,3-2)   <-- this check will fail
 *   }
 *
 *   TEST(test_02)
 *   {
 *       CHECK(1 ,3-2)   <-- this check will pass 
 *       CHECK(2 ,3-2)   <-- this check will fail
 *   }
 *   ...
 *
 * ********************************************************/

struct Test
{
    typedef std::function<void()> testFunc ;

    typedef std::vector<testFunc> TestVec;

    typedef std::map<std::string , TestVec *> TestMap;

    static TestMap & the_map()
    {
        static TestMap themap;
        Test::log.Init("TEST_main");
        return themap;
    }

    static void TRun(const std::string &name ,TestVec * v)
    {
        the_map().emplace(std::make_pair(name,v));
    }

    static void RunAllTest()
    {
        for(auto & a: the_map())
        {
            RunVec(a.first,*a.second);
        }
    }
    static void RunTest(std::string module_name)
    {
        auto it = the_map().find(module_name);
        if(it != the_map().end() )
        {
            RunVec(it->first, *(it->second));
        }
    }
    private:
        static void RunVec(std::string name ,const TestVec & v)
        {
            log<<BGIQD::LOG::lstart()
                <<"---------------- Start test module "<<name<<" ---------------"
                <<BGIQD::LOG::lend();
            for(auto  i = v.begin() ; i!=v.end() ; i++)
            {
                (*i)();
            }
            log<<BGIQD::LOG::lstart()
                <<"---------------- End   test module "<<name<<" ---------------"
                <<BGIQD::LOG::lend();
        }
        static BGIQD::LOG::logger  log;
};


#define TEST_MODULE_INIT(name) \
    static BGIQD::LOG::logger test_logger;\
    static Test::TestVec & get_module()\
    {\
        static Test::TestVec * thevec = nullptr;\
        if(thevec == nullptr)\
        {\
            test_logger.Init("TEST_"#name);\
            thevec = new Test::TestVec();\
            Test::TRun(#name,thevec);\
        }\
        return *thevec;\
    }

#define TEST(name) \
    void name(); \
    namespace{\
        struct test_##name{\
            test_##name(){\
                get_module().push_back([](){\
                    BGIQD::LOG::timer t(test_logger,#name);\
                    name();\
                    });\
            }\
        };\
        static test_##name tmp_##name;\
    }\
    void name()


#define CHECK_STRUCT( expect , value )\
    if( expect != value ) \
        std::cerr<<__FILE__<<":"<<__LINE__<<"struct not match !!"<<std::endl; 

#define CHECK_STRUCT_AND_ONERR( expect , value , action) \
    if( expect != value ) \
    {\
        std::cerr<<__FILE__<<":"<<__LINE__<<"struct not match !!"<<std::endl;\
        action \
    }

#define CHECK_FLOAT( expect , value ) \
    if( ! (expect - value <= 0.000001 || value - expect <= 0.000001) )  \
        std::cerr<<__FILE__<<":"<<__LINE__<<" expect "<<expect<<" but "<<value<<std::endl; 



#define CHECK( expect , value ) \
    if( expect != value ) \
        std::cerr<<__FILE__<<":"<<__LINE__<<" expect "<<expect<<" but "<<value<<std::endl; 

#define CHECK_AND_ONERR( expect , value , action) \
    if( expect != value ) \
    {\
        std::cerr<<__FILE__<<":"<<__LINE__<<" expect "<<expect<<" but "<<value<<std::endl;\
        action \
    }

#endif //__LFR_TEST_TEST_H_
