#include "utils/unittest/Test.h"
#include <unistd.h>
TEST_MODULE_INIT(timetools)

using namespace BGIQD::LOG;

logger unittest_logger;

TEST(timepoint_test)
{
    unittest_logger.Init("timepoint_test");
    timepoint start = timepoint::now();
    sleep(2);
    timepoint end = timepoint::now();
    unittest_logger<<lstart()<<start.to_string()<<lend();
    unittest_logger<<lstart()<<end.to_string()<<lend();
}
TEST(time_period)
{
    unittest_logger.Init("time_period");
    timepoint start = timepoint::now();
    sleep(2);
    timepoint end = timepoint::now();
    unittest_logger<<lstart()<<(end-start).to_string()<<lend();
}

TEST(the_timer)
{
    unittest_logger.Init("the_timer");
    timer(unittest_logger,"test_the_time");
    sleep(2);
}
