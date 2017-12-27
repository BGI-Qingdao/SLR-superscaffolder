#include "Test.h"
#include "argsparser.h"

TEST_MODULE_INIT(Args)

using namespace BGIQD;
using namespace ARGS;

TEST(simulator_args)
{
    auto fake_main=[&](int argc , char ** argv)
    {
        START_PARSE_ARGS

        DEFINE_ARG(bool , open , 'o' )
        DEFINE_ARG(int, k , 'k' )
        DEFINE_ARG(std::string, p , 'p' )

        END_PARSE_ARGS

        CHECK(open.d.b,true);
        CHECK(k.d.i,10);
        CHECK(*(p.d.s) , "/ab/c");
    };
    char a1[] = "main";
    char a2[] = "-o";
    char a3[] = "-k";
    char a4[] = "10";
    char a5[] = "-p";
    char a6[] = "/ab/c";

    char *aa1[100];
    aa1[0]=a1;
    aa1[1]=a2;
    aa1[2]=a3;
    aa1[3]=a4;
    aa1[4]=a5;
    aa1[5]=a6;
    fake_main(6,aa1);
}
