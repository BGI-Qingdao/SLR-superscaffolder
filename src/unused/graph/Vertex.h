#ifndef __GRAPH_VERTEX_H__
#define __GRAPH_VERTEX_H__

#include <vector>

namespace BGIQD{
namespace Graph{

    template<class E>
        struct VertexWithEdge
        {
            typedef E EdgeIndex ;
            std::vector<EdgeIndex> nexts;
        };

}
}

#endif //__GRAPH_VERTEX_H__
