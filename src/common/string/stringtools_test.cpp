#include "Test.h"
#include "stringtools.h"

using namespace BGIQD::STRING;
using namespace BGIQD::LOG;

TEST_MODULE_INIT(String)

TEST(stringsplit)
{
    auto print_vec=[](const std::vector<std::string> &d)
    {
        (test_logger)<<lstart();
        for( size_t i = 0 ; i < d.size() ; i++ ) 
        {
            (test_logger)<<d[i];
        }
        (test_logger)<<lend();
    };

    std::string test = "1234|234";
    std::vector<std::string> expect1;
    expect1.push_back("1234");
    expect1.push_back("234");
    CHECK_STRUCT_AND_ONERR(expect1,split(test,"|"),print_vec(split(test,"|")););
    
    std::string test2 = "1234";
    std::vector<std::string> expect2;
    expect2.push_back("1234");
    CHECK_STRUCT_AND_ONERR(expect2,split(test2,"|"),print_vec(split(test2,"|")););


    std::string test3 = "1234|rs234";
    std::vector<std::string> expect3;
    expect3.push_back("1234");
    expect3.push_back("234");
    CHECK_STRUCT_AND_ONERR(expect3,split(test3,"|rs"),print_vec(split(test3,"|rs")););
}

TEST(stringtrim)
{
    CHECK("ltrim\t\n \r",ltrim("\t\n \rltrim\t\n \r"));
    CHECK("\t\n \rrtrim",rtrim("\t\n \rrtrim\t\n \r"));
    CHECK("ltrim",trim("\t\n \rltrim\t\n \r"));
}
