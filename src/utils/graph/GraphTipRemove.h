#ifndef __ALGORITHM_GRAPH_TIP_REMOVE_H__
#define __ALGORITHM_GRAPH_TIP_REMOVE_H__

#include "utils/graph/GraphBasic.h"
#include <vector>
#include <cassert>
#include <functional>

namespace BGIQD{
    namespace GRAPH{

        /////////////////////////////////////////////////////////////
        //
        // Brief :
        //          Remove tips from a graph iteraterly.
        //
        //          A tip start from tip node and end with junction node.
        //
        //          notice : 
        //              1. it detect all tip .
        //              2. delete tips that shorter than tip_max_num.
        //
        /////////////////////////////////////////////////////////////

        template < class TGraph , class TNode = typename TGraph::Node >
            struct TipRemoveHelper
            {
                typedef TGraph  Graph;
                typedef TNode   Node ;

                typedef typename Graph::NodeId   NodeId ;
                typedef typename Graph::EdgeId   EdgeId ;
                typedef typename Graph::Edge     Edge ;

                typedef std::vector<NodeId> tip ;

                bool debuger ;

                struct TipRemoveResult
                {
                    int base_node_num ;
                    int tip_num ;
                    int tip_node_num ;
                    int base_left_node_num ;
                    std::vector<unsigned int> tip_contigs;

                    void Init()
                    {
                        base_node_num = 0;
                        tip_num = 0;
                        tip_node_num = 0;
                        base_left_node_num = 0;
                        tip_contigs.clear();
                    }
                };

                int tip_max_num;

                TipRemoveResult result ;

                void Init( int max ,bool dd = false )
                {
                    tip_max_num = max;
                    debuger = dd ;
                }

                std::vector<tip>  TipDetect ( const Graph & base )
                {
                    std::vector<tip> ret ;
                    std::set<NodeId> used;
                    for( auto & pair : base.nodes )
                    {
                        auto & x = pair.second ;
                        if( used.find( x.id ) != used.end() )
                        {
                            continue ;
                        }
                        if ( x.EdgeNum() == 1 )
                        {
                            tip tmp_tip ;
                            NodeId curr_id = x.id ;
                            tmp_tip.push_back(x.id);
                            const auto & node = base.GetNode(curr_id) ;
                            typename Node::NodeEdgeIdIterator begin , end;
                            std::tie(begin,end) = node.GetEdges();
                            auto edge_id = *(begin);
                            const auto & edge = base.GetEdge(edge_id);
                            auto next = edge.OppoNode(curr_id) ;
                            auto prev = curr_id ;
                            while( tmp_tip.size() < tip_max_num )
                            {
                                const auto & node = base.GetNode(next) ;
                                if ( node.EdgeNum() == 0 )
                                {
                                    // never should happened .
                                    assert(0);
                                    break ;
                                }
                                else if ( node.EdgeNum() == 1 )
                                {
                                    // find a tip end by tip,so that 
                                    // this is
                                    //      NOT A TIP 
                                    // but a short linear graph.
                                    typename Node::NodeEdgeIdIterator begin , end;
                                    std::tie(begin,end) = node.GetEdges();
                                    auto edge_id = *(begin);
                                    const auto & edge = base.GetEdge(edge_id);
                                    auto next1 = edge.OppoNode(next) ;
                                    assert(next1 == prev );
                                    break ;
                                }
                                else if ( node.EdgeNum() == 2 )
                                {
                                    // tip continue by linear node
                                    tmp_tip.push_back(next);
                                    bool found = false ;
                                    typename Node::NodeEdgeIdIterator begin , end;
                                    std::tie(begin,end) = node.GetEdges();
                                    for( auto x = begin ; x != end ; x++ )
                                    {
                                        const auto & edge = base.GetEdge(*x);
                                        auto next1 = edge.OppoNode(next) ;
                                        if( next1 != prev ) 
                                        {
                                            prev = next ;
                                            next = next1;
                                            found = true ;
                                            break ;
                                        }
                                    }
                                    assert( found );
                                }
                                else
                                {
                                    // tip end by junction node. 
                                    ret.push_back(tmp_tip);
                                    for( auto x : tmp_tip)
                                        used.insert(x);
                                    break;
                                }
                            }
                        }
                        else
                        {
                            ;
                        }
                    }
                    return ret ;
                }

                void TipRemove( Graph & base , const tip & t)
                {
                    for(auto x : t )
                    {
                        if( debuger )
                            std::cerr<<"    del node : "<<x<<std::endl;
                        base.RemoveNode(x);
                    }
                }

                //////////////////////////////////////////////////////
                //
                // Brief :
                //          cyclic remove tips from a graph.
                //
                //          loop :
                //              1. it detect all tip .
                //              2. delete tips that shorter than tip_max_num
                //
                /////////////////////////////////////////////////////
                TipRemoveResult DeepTipRemove( Graph & base )
                {
                    TipRemoveResult ret ;
                    ret.Init() ;
                    ret.base_node_num = base.nodes.size();
                    do
                    {
                        if( debuger )
                            std::cerr<<"    tip round ..."<<std::endl;
                        std::vector<tip> tips ;
                        tips = TipDetect( base );
                        if( tips.empty() )
                        {
                            break;
                        }
                        for( const auto & x : tips )
                        {
                            ret.tip_num ++ ;
                            ret.tip_node_num += x.size() ;
                            TipRemove( base , x);
                            ret.tip_contigs.insert(ret.tip_contigs.begin(), x.begin(),x.end() );
                        }

                    }while( true );
                    ret.base_left_node_num = ret.base_node_num - ret.tip_node_num ;
                    return ret ;
                }
            };

    }
}


#endif //__ALGORITHM_GRAPH_TIP_REMOVE_H__
