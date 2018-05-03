#include "algorithm/graph/DepthSearch.h"
#include "algorithm/graph/Graph.h"
#include "soap2/contigGraph.h"
#include "common/test/Test.h"

TEST_MODULE_INIT(GraphDepthSearch)

//struct TestGraph
//{
//    BGIQD::SOAP2::GraphEA graph_ea ;
//    void Init() 
//    {
//        // init graph
//    }
//};

struct TestGraph
{
    struct Node
    {
        char id ;
        int edge_id ;
    };

    struct Edge
    {
        char from ; 
        char to ; 
        int next ;
    };

    std::map< char , Node> nodes ;
    std::map< int , Edge> edges ;

    /**************************************************
     * Test graph :
     *
     * digraph{
     *  'u'->'v'
     *  'u'->'x'
     *  'x'->'v'
     *  'v'->'y'
     *  'y'->'x'
     *  'w'->'y'
     *  'w'->'z'
     *  'z'->'z'
     * }
     *
     * draw the graph by DOT if you want.
     *************************************************/

    static TestGraph GetTestGraph()
    {

        TestGraph test;

        test.add_edge(1,'u','v',2);
        test.add_edge(2,'u','x',-1);
        test.add_edge(3,'x','v',-1);
        test.add_edge(4,'v','y',-1);
        test.add_edge(5,'y','x',-1);
        test.add_edge(6,'w','y',7);
        test.add_edge(7,'w','z',-1);
        test.add_edge(8,'z','z',-1);

        test.add_node('u',1);
        test.add_node('x',3);
        test.add_node('v',4);
        test.add_node('y',5);
        test.add_node('w',6);
        test.add_node('z',8);

        return test ;
    }

    private :

    void add_edge( int id , char from , char to , int next )
    {
        auto & edge3 = edges[id] ;
        edge3.from = from;
        edge3.to= to;
        edge3.next = next ;
    }

    void add_node( char id , int edge_id )
    {
        auto & node = nodes[id];
        node.id = id ;
        node.edge_id  = edge_id;
    }

};

template<>
struct BGIQD::GRAPH::GraphAccess<TestGraph , char , int >
{
    typedef char                         GraphNodeId ;
    typedef int                          GraphEdgeId ;

    typedef GraphNodeBase<char ,int>     Node;
    typedef GraphEdgeBase<char ,int>     Edge;

    Node & AccessNode(char i)
    {
        auto itr = nodes.find(i); 
        if( itr == nodes.end() )
        {
            auto & n = nodes[i] ;
            n.id = (*base).nodes.at(i).id ;
            n.edge_id = (*base).nodes.at(i).edge_id ;
            return n;
        }
        return itr->second ;
    }

    Edge & AccessEdge(int i)
    {
        auto itr = edges.find(i); 
        if( itr == edges.end() )
        {
            auto & n = edges[i] ;
            auto & base_edge = base->edges.at(i);
            n.id = i ;
            n.from = base_edge.from;
            n.to   = base_edge.to;
            n.next = base_edge.next;
            return n;
        }
        return itr->second ;
    }

    TestGraph * base ;
    std::map< char, Node> nodes ;
    std::map< int , Edge> edges ;
};

struct traits_1 {} ;
struct traits_2 {} ;

typedef BGIQD::GRAPH::GraphAccess<TestGraph,char,int> ACCESS;
typedef BGIQD::GRAPH::EdgeIterator<ACCESS> EdgeItr;

template<>
struct BGIQD::GRAPH::DepthSearchPathEndHelper<ACCESS, traits_1>
{
    typedef typename ACCESS::Node            Node;
    typedef typename ACCESS::Edge            Edge;
    typedef typename ACCESS::GraphNodeId     NodeId;
    typedef typename ACCESS::GraphEdgeId     EdgeId;

    typedef traits_1                            traisId;
    void AddNode(const Node & ) {} 
    void AddEdge(const Edge & ) {
    }

    void PopEdge() {}
    void PopNode() {}

    bool IsEnd() const {
        //
    }
    
    void Init( ACCESS & acc) { 
        accesser = & acc ;
    }
    private:
    ACCESS * accesser;
    Node curr 
    //std::stack<Node> nodepath;
    //std::stack<Edge> edgepath;
};


TEST(GraphAccessNode)
{
    auto test = TestGraph::GetTestGraph() ;

    BGIQD::GRAPH::GraphAccess<TestGraph,char,int> acc ;
    acc.base = &test;

    CHECK('u',acc.AccessNode('u').id);
    CHECK(1,acc.AccessNode('u').edge_id);

    CHECK('v',acc.AccessNode('v').id);
    CHECK(4,acc.AccessNode('v').edge_id);

    CHECK('x',acc.AccessNode('x').id);
    CHECK(3,acc.AccessNode('x').edge_id);

    CHECK('y',acc.AccessNode('y').id);
    CHECK(5,acc.AccessNode('y').edge_id);

    CHECK('w',acc.AccessNode('w').id);
    CHECK(6,acc.AccessNode('w').edge_id);

    CHECK('z',acc.AccessNode('z').id);
    CHECK(8,acc.AccessNode('z').edge_id);
}

TEST(GraphAccessEdge)
{
    auto test = TestGraph::GetTestGraph() ;

    BGIQD::GRAPH::GraphAccess<TestGraph,char,int> acc ;
    acc.base = &test;

    CHECK(1  ,acc.AccessEdge(1).id);
    CHECK('u',acc.AccessEdge(1).from);
    CHECK('v',acc.AccessEdge(1).to);
    CHECK(2  ,acc.AccessEdge(1).next);

    CHECK(2  ,acc.AccessEdge(2).id);
    CHECK('u',acc.AccessEdge(2).from);
    CHECK('x',acc.AccessEdge(2).to);
    CHECK(-1 ,acc.AccessEdge(2).next);

    CHECK(5  ,acc.AccessEdge(5).id);
    CHECK('y',acc.AccessEdge(5).from);
    CHECK('x',acc.AccessEdge(5).to);
    CHECK(-1 ,acc.AccessEdge(5).next);

    CHECK(8  ,acc.AccessEdge(8).id);
    CHECK('z',acc.AccessEdge(8).from);
    CHECK('z',acc.AccessEdge(8).to);
    CHECK(-1 ,acc.AccessEdge(8).next);
}

TEST(EdgeItr_test)
{
    auto test = TestGraph::GetTestGraph() ;
    ACCESS acc ;
    acc.base = &test;

    EdgeItr itr1(acc.AccessEdge(1),acc);
    CHECK( 1 , itr1->id );
    CHECK('u', itr1->from);
    CHECK('v', itr1->to);
    CHECK(2  , itr1->next);
    ++ itr1 ;

    CHECK( 2 , itr1->id );
    CHECK('u', itr1->from);
    CHECK('x', itr1->to);
    CHECK(-1  , itr1->next);
    ++ itr1 ;

    CHECK( true , (EdgeItr::end() == itr1 ));

    EdgeItr itr2(acc.AccessEdge(8),acc);
    CHECK( 8 , itr2->id );
    CHECK('z', itr2->from);
    CHECK('z', itr2->to);
    CHECK(-1 , itr2->next);
    ++ itr2 ;

    CHECK( true , (EdgeItr::end() == itr2 ));
}


TEST(GraphDepthSearch_test)
{

}
