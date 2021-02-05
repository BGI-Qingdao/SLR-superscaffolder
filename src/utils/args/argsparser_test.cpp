#include "utils/unittest/Test.h"
#include "argsparser.h"

TEST_MODULE_INIT(Args)

using namespace BGIQD;
using namespace ARGS;

TEST(simulator_args)
{
    auto fake_main=[&](int argc , char ** argv)
    {
        START_PARSE_ARGS

        DEFINE_ARG_REQUIRED(bool , open , "require open" );
        DEFINE_ARG_REQUIRED(int, k , "require k" );
        DEFINE_ARG_REQUIRED(float , f , "require f" );
        DEFINE_ARG_OPTIONAL(std::string, p , "optional p", "default p" );
        DEFINE_ARG_OPTIONAL(std::string, oo , "optional oo", "ooo" );
        DEFINE_ARG_OPTIONAL(bool, tb, "optional tb", "" );
        DEFINE_ARG_OPTIONAL(float, tf, "optional tf", "1.0f" );
        DEFINE_ARG_OPTIONAL(int, ti, "optional tf", "10" );

        END_PARSE_ARGS

        // check assigned parameters:
        CHECK(true,open.to_bool());
        CHECK(10,k.to_int());
        CHECK("/ab/c",p.to_string());
        CHECK(5.0,f.to_float());
        // check default parameters:
        CHECK("ooo" ,oo.to_string() );
        CHECK(false, tb.to_bool());
        CHECK(1.0f, tf.to_float());
        CHECK(10,ti.to_int());
        return 1;
    };
    char a1[] = "main";
    char a2[] = "-open";
    char a3[] = "--k";
    char a4[] = "10";
    char a5[] = "--p";
    char a6[] = "/ab/c";
    char a7[] = "--f";
    char a8[] = "5.0";

    char *aa1[10];
    aa1[0]=a1;
    aa1[1]=a2;
    aa1[2]=a3;
    aa1[3]=a4;
    aa1[4]=a5;
    aa1[5]=a6;
    aa1[6]=a7;
    aa1[7]=a8;

    fake_main(8,aa1);

    char *aa2[1];
    aa2[0]=a1;
    // will print usage :
    fake_main(1,aa2);
}
