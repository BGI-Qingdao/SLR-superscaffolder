#ifndef __LFR_TEST_TEST_H__
#define __LFR_TEST_TEST_H__

#include "Check.h"
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <functional>
#include "log.h"
#include "logfilter.h"
struct Test
{
    typedef std::function<void()> testFunc ;

    typedef std::vector<testFunc> TestVec;

    typedef std::map<std::string , TestVec *> TestMap;

    static TestMap & the_map()
    {
        static TestMap themap;
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
            std::cout<<"---------------- Start test module "<<name<<" ---------------"<<std::endl;
            for(auto  i = v.begin() ; i!=v.end() ; i++)
            {
                (*i)();
            }
            std::cout<<"---------------- End   test module "<<name<<" ---------------"<<std::endl<<std::endl;
        }
};


#define TEST_MODULE_INIT(name) \
    static BGIQD::LOG::logger * test_logger;\
    static Test::TestVec & get_module()\
    {\
        static Test::TestVec * thevec = nullptr;\
        if(thevec == nullptr)\
        {\
            test_logger = BGIQD::LOG::logfilter::singleton().get("TEST"#name,BGIQD::LOG::loglevel::DEBUG);\
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


#endif //__LFR_TEST_TEST_H_
