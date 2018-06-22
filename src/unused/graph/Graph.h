#ifndef __GRAPH_GRAPH_H__
#define __GRAPH_GRAPH_H__

#include <vector>
#include <map>

#include "graph/Edge.h"
#include "graph/Vertex.h"

namespace BGIQD{
namespace Graph{

    template<class I , class E >
        class DiGraph_VE
        {
            public:
                typedef I VertexIndex;
                typedef E EdgeIndex;
                typedef VertexWithEdge<E> G_Vertex;
                typedef Edge<I>  G_Edge;

                std::vector<G_Vertex> vertexs;
                std::vector<G_Edge> edges;

            public:
                std::map<std::pair<I,I> ,E> connections; 
            private:
                EdgeIndex nextE;
                VertexIndex nextV;


            public:

                EdgeIndex AddEdge( I from , I to )
                {
                    EdgeIndex curr = edges.size();
                    edges.emplace_back(from,to);

                    vertexs[from].nexts.push_back(curr);

                    connections[std::make_pair(from,to)] =curr;

                    return curr ;
                }

                bool IsConnected(I i1 , I i2 , E & e)
                {
                    auto itr = connections.find(std::make_pair(i1 ,i2));
                    if( itr != connections.end() )
                    {
                        e = itr->second;
                        return true;
                    }
                    return false;
                }

                VertexIndex From( EdgeIndex e)
                {
                    return edges[e-1].from;
                }
                VertexIndex To( EdgeIndex e)
                {
                    return edges[e-1].to;
                }
        };


        struct ContigVE
        {
            
        };

}
}

#endif //__GRAPH_GRAPH_H__
