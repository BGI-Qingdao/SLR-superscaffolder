#ifndef __GRAPH_EDGE_H__
#define __GRAPH_EDGE_H__
namespace BGIQD{
namespace Graph{

    template<class V>
        class Edge{
            public:
                typedef V VertexIndex;

                VertexIndex from;
                VertexIndex to;
        };
}
}

#endif //__GRAPH_EDGE_H_
