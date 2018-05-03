#ifndef __ALGORITHM_GRAPH_GRAPH_H__
#define __ALGORITHM_GRAPH_GRAPH_H__ 

#include <iterator>

namespace BGIQD {
    namespace GRAPH {

        template<class NodeId ,  class EdgeId >
            struct GraphEdgeBase
            {
                typedef NodeId EdgeNodeId;
                typedef EdgeId EdgeEdgeId;

                EdgeId      id ;
                EdgeId      next;

                NodeId      from ;
                NodeId      to ;

                static const EdgeId invalid = -1 ;
            } ;



        template<class NodeId ,  class EdgeId >
            struct GraphNodeBase
            {
                typedef NodeId NodeNodeId;
                typedef EdgeId NodeEdgeId;

                NodeId                          id ;
                EdgeId                          edge_id;
            };


        // Must be specialied before use
        template<class BaseGraph , class NodeId ,  class EdgeId>
            struct GraphAccess
            {
                typedef NodeId                       GraphNodeId ;
                typedef EdgeId                       GraphEdgeId ;

                typedef GraphNodeBase<NodeId,EdgeId> Node;
                typedef GraphEdgeBase<NodeId,EdgeId> Edge;

                Node & AccessNode(NodeId);
                Edge & AccessEdge(EdgeId);

            };

        template<class GraphAccess>
            struct EdgeIterator : public std::iterator<std::forward_iterator_tag,int>
        {
            public:

                typedef typename GraphAccess::GraphEdgeId Id;
                typedef typename GraphAccess::Edge        Edge;

                EdgeIterator() : curr(NULL) ,accessor( NULL ) { }

                EdgeIterator(const Edge & e , GraphAccess & acc )
                {
                    curr = &e;
                    accessor = &acc ;
                }

                EdgeIterator( const EdgeIterator & ei )
                {
                    curr = ei.curr ;
                    accessor = ei.accessor ;
                }

                EdgeIterator & operator = ( const EdgeIterator & ei )
                {
                    if( &ei != this )
                    {
                        curr = ei.curr ;
                        accessor = ei.accessor ;
                    }
                    return *this;
                }

                bool operator == ( const EdgeIterator & ei )const
                {
                    return curr == ei.curr ;
                }

                bool operator != ( const EdgeIterator & ei )const
                {
                    return curr != ei.curr ;
                }

                const Edge & operator*() const  { return *curr ; }

                const Edge * operator->() const  { return curr ; }

                EdgeIterator & operator ++() {
                    if( curr != NULL && accessor != NULL )
                    {
                        Id next = curr->next ;
                        if( next != Edge::invalid )
                            curr = &(accessor->AccessEdge(next));
                        else
                            curr = NULL ;
                    }
                    else
                        curr = NULL ;
                    return *this ;
                }

                static EdgeIterator & end() 
                {
                    static EdgeIterator end;
                    return end ;
                }
            private:

                const Edge * curr ;
                GraphAccess * accessor ;
        };

    } // namespace GRAPH
} // namespace BGIQD

#endif //__ALGORITHM_GRAPH_GRAPH_H__
