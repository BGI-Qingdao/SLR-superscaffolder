#include "Test.h"
#include "argsparser.h"

TEST_MODULE_INIT(Args)

using namespace BGIQD;
using namespace ARGS;

TEST(simulator_args)
{
    auto fake_main=[](int argc , char ** argv)
    {
        START_PARSE_ARGS

        DEFINE_ARG(bool , open , 'o' )

        END_PARSE_ARGS

    };
}
