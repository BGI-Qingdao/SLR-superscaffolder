#include "utils/unittest/Test.h"
#include "utils/graph/GraphBasic.h"
#include "utils/graph/mst/MinTree.h"

TEST_MODULE_INIT(MinTree)


typedef BGIQD::GRAPH::IGraphNodeBasic<std::string , int > MTestNode ;

typedef BGIQD::GRAPH::IGraphEdgeBasic<std::string , int > MTestEdge;

struct MTEdge : public MTestEdge
{
    int value ;
    std::string ToString() const
    {
        std::ostringstream ost;
        ost<<from<<"\t--\t"<<to<<" [ id = "<<id << " value = "<<value<< " ]";
        return ost.str();
    }
};

struct EAttr
{
    int GetValue(const MTEdge & e ) const
    {
        return e.value ;
    }
};


struct MTestGraphBasic : public BGIQD::GRAPH::Graph<MTestNode , MTEdge>
{
    typedef BGIQD::GRAPH::Graph<MTestNode , MTEdge> Basic;

    void AddEdgeValue( const std::string & from , const std::string & to, int value )
    {
        MTEdge tmp ;
        tmp.from = from ;
        tmp.to = to ;
        tmp.value = value ;
        Basic::AddEdge(tmp);
    }

    void AddNode(const MTestNode & n){
        Basic::AddNode(n);
    }
    void AddNode(const std::string & str){
        MTestNode tmp;
        tmp.id = str;
        AddNode(tmp);
    }

    static MTestGraphBasic TestData() {
        MTestGraphBasic test ;
        test.AddNode("a");
        test.AddNode("b");
        test.AddNode("c");
        test.AddNode("d");
        test.AddNode("e");
        test.AddNode("f");
        test.AddNode("g");
        test.AddNode("h");
        test.AddNode("i");
        test.AddNode("i");
        test.AddEdgeValue( "a" , "b",4);
        test.AddEdgeValue( "a" , "h",8);
        test.AddEdgeValue( "b" , "c",8);
        test.AddEdgeValue( "b" , "h",11);
        test.AddEdgeValue( "c" , "d",7);
        test.AddEdgeValue( "c" , "f",4);
        test.AddEdgeValue( "h" , "i" ,7);
        test.AddEdgeValue( "h" , "g" ,1);
        test.AddEdgeValue( "i" , "c" ,2);
        test.AddEdgeValue( "i" , "g" ,6);
        test.AddEdgeValue( "g" , "f" ,2);
        test.AddEdgeValue( "f" , "d" ,17);
        test.AddEdgeValue( "d" , "e" ,9);
        test.AddEdgeValue( "f" , "e" ,10);
        return test;
    }

};

typedef BGIQD::GRAPH::MinTreeHelper<MTestGraphBasic, int ,EAttr> MTHelper;

/********************************************************************
*
* The initial graph       ---->>>    The mst
*   +--------------+                  +--------------+
*   |              |                  |              |
*   |           +4-f-17-+             |           +4-f
*   |           |  |    |             |           |
*   2           |  10   |             2           | 
*   |           |  |    |             |           |
*   |           |  e-+  |             |           |  e-+
*   |           |    9  |             |           |    9
*   |           |    |  |             |           |    |
*   |   a-4-b-8-c-7--d--+             |   a-4-b   c-7--d
*   |   |   |   |                     |   |       |
*   |   8   11  |                     |   8       |
*   |   |   |   |                     |   |       |
*   g-1-h---+   |                     g-1-h       |
*   |   |       |                                 |
*   |   7---i-2-+                             i-2-+
*   |       |
*   +--6----+
*
*
* *****************************************************************/
TEST(MinTree)
{
    auto test = MTestGraphBasic::TestData();
    MTHelper mtHelper;
    EAttr attr;
    auto m = mtHelper.MinTree(test,attr);
    m.PrintAsDOT(std::cout);
    CHECK(true, m.CheckEdge("a","b"));
    CHECK(true, m.CheckEdge("a","h"));
    CHECK(true, m.CheckEdge("h","g"));
    CHECK(true, m.CheckEdge("g","f"));
    CHECK(true, m.CheckEdge("c","f"));
    CHECK(true, m.CheckEdge("i","c"));
    CHECK(true, m.CheckEdge("d","e"));
    CHECK(true, m.CheckEdge("c","d"));

    CHECK(false, m.CheckEdge("b","c"));
    CHECK(false, m.CheckEdge("b","h"));
    CHECK(false, m.CheckEdge("h","i"));
    CHECK(false, m.CheckEdge("i","g"));
    CHECK(false, m.CheckEdge("f","d"));
    CHECK(false, m.CheckEdge("f","e"));
}
