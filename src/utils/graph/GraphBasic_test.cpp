#include "utils/unittest/Test.h"
#include "utils/graph/GraphBasic.h"

TEST_MODULE_INIT(GraphBasicTest)

typedef BGIQD::GRAPH::IGraphNodeBasic<std::string , int > TestNode ;

TEST(Test_Node)
{
    TestNode t1;
    t1.AddEdge(1);
    t1.AddEdge(2);
    t1.AddEdge(3);
    CHECK(3,t1.EdgeNum());
    CHECK(true,t1.HasEdge(1));
    t1.RemoveEdge(1);
    CHECK(false,t1.HasEdge(1));
    CHECK(2,t1.EdgeNum());
    CHECK(true,t1.HasEdge(2));
    t1.CleanEdges();
    CHECK(false,t1.HasEdge(1));
    CHECK(false,t1.HasEdge(1));
    CHECK(0,t1.EdgeNum());
    CHECK(false,t1.HasEdge(2));
}

TEST(Test_edge_iterator)
{
    TestNode t1;
    int i ;
    for( i = 1 ; i < 100 ; i ++ )
        t1.AddEdge(i);
    typedef TestNode::NodeEdgeIdIterator Itr;
    Itr begin,end;
    std::tie(begin,end) = t1.GetEdges();
    i = 1 ;
    for( Itr x = begin ; x != end ; x++ , i++){
        CHECK(i,*x)
    }
}

typedef BGIQD::GRAPH::IGraphEdgeBasic<std::string , int > TestEdge;

TEST( Test_Edge )
{
    TestEdge t1;
    t1.from = "a";t1.to="b";
    TestEdge t2;

    t2.from = "a";t2.to="b";
    bool eq=(t1==t2);
    bool neq=(t1!=t2);
    CHECK(true, eq);
    CHECK(false, neq);

    t2.to = "c";
    eq=(t1==t2);
    neq=(t1!=t2);
    CHECK(false, eq);
    CHECK(true, neq);

    CHECK("a",t1.OppoNode("b"));
    CHECK("b",t1.OppoNode("a"));
    CHECK("a",t2.OppoNode("c"));
    CHECK("c",t2.OppoNode("a"));
    t1.InvalidMe();
    CHECK(false, t1.IsValid());
    CHECK(true, t2.IsValid());
}

typedef BGIQD::GRAPH::IDigraphEdgeBase<std::string , int > TestDiEdge;

TEST( Test_DiEdge )
{
    TestDiEdge t1;
    t1.from = "a";t1.to="b";
    TestDiEdge t2;
    t2.from = "a";t2.to="b";

    bool eq=(t1==t2);
    bool neq=(t1!=t2);
    CHECK(true, eq);
    CHECK(false, neq);

    t2.to = "c";
    eq=(t1==t2);
    neq=(t1!=t2);
    CHECK(false, eq);
    CHECK(true, neq);

    t2.from = "b";t2.to="a";
    eq=(t1==t2);
    neq=(t1!=t2);
    CHECK(false, eq);
    CHECK(true, neq);

    CHECK("a",t1.OppoNode("b"));
    CHECK("b",t1.OppoNode("a"));
    CHECK("b",t2.OppoNode("a"));
    CHECK("a",t2.OppoNode("b"));
    t1.InvalidMe();
    CHECK(false, t1.IsValid());
    CHECK(true, t2.IsValid());
}

struct TestGraphBasic : public BGIQD::GRAPH::Graph<TestNode , TestEdge>
{
    typedef BGIQD::GRAPH::Graph<TestNode , TestEdge> Basic;

    void AddNode(const std::string & name) {
        TestNode n1;
        n1.id = name ;
        Basic::AddNode(n1);
    }

    void AddEdge( const std::string &N1, const std::string &N2) {
        TestEdge tmp;
        tmp.from = N1 ;
        tmp.to = N2;
        Basic::AddEdge(tmp);
    }

    static TestGraphBasic TestData() {
        TestGraphBasic test ;
        test.AddNode("NodeA");
        test.AddNode("NodeB");
        test.AddNode("NodeC");
        test.AddNode("NodeD");
        test.AddEdge( "NodeA" , "NodeB");
        test.AddEdge( "NodeA" , "NodeC");
        test.AddEdge( "NodeB" , "NodeC");
        test.AddEdge( "NodeB" , "NodeD");
        test.AddEdge( "NodeC" , "NodeD");
        return test;
    }
};

struct TestDigraphBasic : public BGIQD::GRAPH::Digraph<TestNode , TestDiEdge>
{
    typedef BGIQD::GRAPH::Digraph<TestNode , TestDiEdge> Basic;

    void AddNode(const std::string & name) {
        TestNode n1;
        n1.id = name ;
        Basic::AddNode(n1);
    }

    void AddEdge( const std::string &N1, const std::string &N2) {
        TestDiEdge tmp;
        tmp.from = N1 ;
        tmp.to = N2;
        Basic::AddEdge(tmp);
    }

    static TestDigraphBasic TestData() {
        TestDigraphBasic test ;
        test.AddNode("NodeA");
        test.AddNode("NodeB");
        test.AddNode("NodeC");
        test.AddNode("NodeD");
        test.AddEdge( "NodeA" , "NodeB");
        test.AddEdge( "NodeA" , "NodeC");
        test.AddEdge( "NodeB" , "NodeC");
        test.AddEdge( "NodeB" , "NodeD");
        test.AddEdge( "NodeC" , "NodeD");
        return test;
    }
};


TEST(ListGraph)
{
    auto  test = TestGraphBasic::TestData();
    test.PrintAsDOT(std::cout);
    CHECK( true ,test.CheckEdge("NodeA","NodeB") ) 
    CHECK( true ,test.CheckEdge("NodeB","NodeA") ) 
    CHECK( true , test.GetNode("NodeA").HasEdge(0));
    CHECK( true , test.GetNode("NodeB").HasEdge(0));
}

TEST(ListDigraph)
{
    auto  test = TestDigraphBasic::TestData();
    test.PrintAsDOT(std::cout);
    CHECK( true ,test.CheckEdge("NodeA","NodeB") ) 
    CHECK( false ,test.CheckEdge("NodeB","NodeA") ) 
    CHECK( true , test.GetNode("NodeA").HasEdge(0));
    CHECK( false, test.GetNode("NodeB").HasEdge(0));
}

