#include "utils/unittest/Test.h"
#include "utils/graph/DisJointSet.h"

TEST_MODULE_INIT(DisJoinSetTest)

using namespace BGIQD::GRAPH;
typedef DisJoin_Set<char>  DS;
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
    test.AddConnect('A','B');
    test.AddConnect('B','C');
    test.AddConnect('C','D');
    test.AddConnect('A','F');
    test.AddConnect('B','D');
    test.AddConnect('X','Y');
    test.AddConnect('X','Z');
    test.AddConnect('O','X');
    test.AddConnect('Z','F');
    test.AddConnect('M','N');

    CHECK('X',test.GetGroup('A'));
    CHECK('X',test.GetGroup('B'));
    CHECK('X',test.GetGroup('C'));
    CHECK('X',test.GetGroup('D'));
    CHECK('X',test.GetGroup('F'));
    CHECK('X',test.GetGroup('X'));
    CHECK('X',test.GetGroup('Y'));
    CHECK('X',test.GetGroup('O'));
    CHECK('X',test.GetGroup('Z'));

    CHECK('M',test.GetGroup('M'));
    CHECK('M',test.GetGroup('N'));
}
