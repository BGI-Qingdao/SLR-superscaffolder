#include "common/test/Test.h"
#include "algorithm/graph/GraphTipRemove.h"
#include "algorithm/graph/GraphBasic.h"

TEST_MODULE_INIT(GraphTipRemove)

typedef BGIQD::GRAPH::IGraphNodeBasic<std::string , int > TipTestNode ;

typedef BGIQD::GRAPH::IGraphEdgeBasic<std::string , int > TipTestEdge;

struct TipTestGraphBasic : public BGIQD::GRAPH::ListGraph<TipTestNode , TipTestEdge>
{
    typedef BGIQD::GRAPH::ListGraph<TipTestNode, TipTestEdge> Basic;

    void AddEdgeValue( const std::string & from , const std::string & to)
    {
        TipTestEdge tmp ;
        tmp.from = from ;
        tmp.to = to ;
        Basic::AddEdge(tmp);
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
         *
         * after 2 round of tip( detect only 1 node as tip ) remove :
         *
         * left :
         *       g---h---i--j--k
         *                     |
         *                     +---n---o
         */

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
        return test;
    }

};

typedef BGIQD::GRAPH::TipRemoveHelper<TipTestGraphBasic> TipHelper;

bool tip_checker(const TipHelper::tip & t )
{
    return t.size() < 2 ;
}

TEST(TipRemove_test1)
{

    auto test_graph = TipTestGraphBasic::TestData1();

    TipHelper tester;

    tester.Init(tip_checker,true);
    std::cerr<<"before tip remove"<<std::endl;
    test_graph.PrintAsDOT(std::cerr);
    tester.DeepTipRemove(test_graph);
    std::cerr<<"after tip remove"<<std::endl;
    test_graph.PrintAsDOT(std::cerr);
    std::cerr<<std::endl;
}

