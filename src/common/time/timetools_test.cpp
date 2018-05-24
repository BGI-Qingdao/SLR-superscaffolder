#include "common/test/Test.h"
#include "common/time/timetools.h"
#include "common/test/Check.h"
#include <unistd.h>

TEST_MODULE_INIT(timetools)

using namespace BGIQD::TIME;
using namespace BGIQD::LOG;

TEST(timepoint_test)
{
    timepoint start = timepoint::now();
    sleep(2);
    timepoint end = timepoint::now();
    (test_logger)<<lstart()<<start.to_string()<<lend();
    (test_logger)<<lstart()<<end.to_string()<<lend();
    (test_logger)<<lstart()<<(end-start).to_string()<<lend();
}


