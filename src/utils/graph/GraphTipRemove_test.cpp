#include "utils/unittest/Test.h"
#include "utils/graph/GraphTipRemove.h"
#include "utils/graph/GraphBasic.h"

TEST_MODULE_INIT(GraphTipRemove)

typedef BGIQD::GRAPH::IGraphNodeBasic<std::string , int > TipTestNode ;

typedef BGIQD::GRAPH::IGraphEdgeBasic<std::string , int > TipTestEdge;

struct TipTestGraphBasic : public BGIQD::GRAPH::Graph<TipTestNode , TipTestEdge>
{
    typedef BGIQD::GRAPH::Graph<TipTestNode, TipTestEdge> Basic;

    void AddEdgeValue( const std::string & from , const std::string & to)
    {
        TipTestEdge tmp ;
        tmp.from = from ;
        tmp.to = to ;
        Basic::AddEdge(tmp);
    }
    void AddNode(const std::string & c ){
        TipTestNode tmp;
        tmp.id = c ;
        Basic::AddNode(tmp);
    }
    static TipTestGraphBasic TestData1() {
        TipTestGraphBasic test ;

        /*
         * Init graph :
         *
         * a---+
         *     e--+             +------m
         * b---+  |             |
         *        g---h---i--j--k 
         * c --+  |             |
         *     f--+             +------n--o
         * d --+
         * aa --- bb
         *
         * Round 1 will delete a,b,c,d,m :
         *
         *     e--+
         *        |
         *        g---h---i--j--k 
         *        |             |
         *     f--+             +------n--o
         * aa --- bb
         *
         * Round 2 will delete e,f :
         *
         *        g---h---i--j--k 
         *                      |
         *                      +------n--o
         * aa --- bb
         *
         * END rip remove
         */
        test.AddNode("a");
        test.AddNode("b");
        test.AddNode("c");
        test.AddNode("d");
        test.AddNode("e");
        test.AddNode("f");
        test.AddNode("g");
        test.AddNode("h");
        test.AddNode("i");
        test.AddNode("j");
        test.AddNode("k");
        test.AddNode("m");
        test.AddNode("n");
        test.AddNode("o");
        test.AddNode("aa");
        test.AddNode("bb");
        test.AddEdgeValue( "a" , "e");
        test.AddEdgeValue( "b" , "e");
        test.AddEdgeValue( "c" , "f");
        test.AddEdgeValue( "d" , "f");
        test.AddEdgeValue( "e" , "g");
        test.AddEdgeValue( "f" , "g");
        test.AddEdgeValue( "g" , "h");
        test.AddEdgeValue( "h" , "i");
        test.AddEdgeValue( "i" , "j");
        test.AddEdgeValue( "j" , "k");
        test.AddEdgeValue( "k" , "m");
        test.AddEdgeValue( "k" , "n");
        test.AddEdgeValue( "n" , "o");
        test.AddEdgeValue( "aa" , "bb");
        return test;
    }

};

typedef BGIQD::GRAPH::TipRemoveHelper<TipTestGraphBasic> TipHelper;

TEST(TipRemove_test1)
{

    auto test_graph = TipTestGraphBasic::TestData1();

    TipHelper tester;

    tester.Init(2,true);
    std::cerr<<"before tip remove"<<std::endl;
    test_graph.PrintAsDOT(std::cerr);
    auto ret = tester.DeepTipRemove(test_graph);
    CHECK(16 , ret.base_node_num );
    CHECK( 9 , ret.base_left_node_num );
    CHECK( 7 , ret.tip_num );
    CHECK( 7 , ret.tip_node_num ) ;
    CHECK(false, test_graph.HasNode("a"));
    CHECK(false, test_graph.HasNode("b"));
    CHECK(false, test_graph.HasNode("c"));
    CHECK(false, test_graph.HasNode("d"));
    CHECK(false, test_graph.HasNode("e"));
    CHECK(false, test_graph.HasNode("f"));
    CHECK(false, test_graph.HasNode("m"));
    CHECK(true, test_graph.HasNode("g"));
    CHECK(true, test_graph.HasNode("h"));
    CHECK(true, test_graph.HasNode("i"));
    CHECK(true, test_graph.HasNode("g"));
    CHECK(true, test_graph.HasNode("k"));
    CHECK(true, test_graph.HasNode("n"));
    CHECK(true, test_graph.HasNode("o"));
    CHECK(true, test_graph.HasNode("aa"));
    CHECK(true, test_graph.HasNode("bb"));
    std::cerr<<"after tip remove"<<std::endl;
    test_graph.PrintAsDOT(std::cerr);
    std::cerr<<std::endl;
}
