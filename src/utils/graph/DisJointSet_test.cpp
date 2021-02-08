#include "utils/unittest/Test.h"
#include "utils/graph/DisJointSet.h"

TEST_MODULE_INIT(DisJoinSetTest)

using namespace BGIQD::GRAPH;
typedef DisJoin_Set  DS;
static int A=0;
static int B=1;
static int C=2;
static int D=3;
static int E=4;
static int F=5;
static int X=6;
static int Y=7;
static int Z=8;
static int O=9;
static int M=10;
static int N=11;
/**********************************
 * the data:    | the graph :
 * set1:        |
 *  A-B         |    F-A-B-C-D
 *  B-C         |    |   |___|
 *  C-D         |    |
 *  A-F         |    |
 *  B-D         |    Z-X-Y
 *              |      |
 *  X-Y         |      O
 *  X-Z         |
 *  O-X         |
 *  Z-F         |
 *              |
 * set2:        |
 *  M-N         |   M-N
 *
 * ******************************/
TEST(Test_disjoin)
{
    DS test;
    test.AddConnect(A,B);
    test.AddConnect(B,C);
    test.AddConnect(C,D);
    test.AddConnect(A,F);
    test.AddConnect(B,D);
    test.AddConnect(X,Y);
    test.AddConnect(X,Z);
    test.AddConnect(O,X);
    test.AddConnect(Z,F);
    test.AddConnect(M,N);
    test.GenAllResult();
    CHECK(0,test.GetGroup(A));
    CHECK(0,test.GetGroup(B));
    CHECK(0,test.GetGroup(C));
    CHECK(0,test.GetGroup(D));
    CHECK(0,test.GetGroup(F));
    CHECK(0,test.GetGroup(X));
    CHECK(0,test.GetGroup(Y));
    CHECK(0,test.GetGroup(O));
    CHECK(0,test.GetGroup(Z));

    CHECK(2,test.GetGroup(M));
    CHECK(2,test.GetGroup(N));
}
