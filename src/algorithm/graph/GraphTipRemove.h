#ifndef __ALGORITHM_GRAPH_TIP_REMOVE_H__
#define __ALGORITHM_GRAPH_TIP_REMOVE_H__

#include "algorithm/graph/GraphBasic.h"
#include <vector>
#include <cassert>
#include <functional>

namespace BGIQD{
    namespace GRAPH{


        //////////////////////////////////////////////////////
        //
        // Brief :
        //          Remove tips from a graph.
        //
        //          notice : 
        //              1. it is detect all tip .
        //              2. delete tips that shorter than tip_cut_num
        //
        /////////////////////////////////////////////////////
        template < class TListGraph , class TNode = typename TListGraph::Node >
            struct TipRemoveHelper
            {
                typedef TListGraph  Graph;
                typedef TNode       Node ;

                typedef typename Graph::NodeId   NodeId ;
                typedef typename Graph::EdgeId   EdgeId ;

                typedef std::vector<NodeId> tip;

                // return true if tip is still short enough.
                //        false otherwise .
                typedef std::function<bool(const tip & )> TipChecker;

                struct TipRemoveResult
                {
                    int base_node_num ;
                    int tip_num ;
                    int tip_node_num ;
                    int base_left_node_num ;

                    void Init()
                    {
                        base_node_num = 0;
                        tip_num = 0;
                        tip_node_num = 0;
                        base_left_node_num = 0;
                    }
                };

                TipChecker  checker;

                TipRemoveResult result ;

                void Init( TipChecker c )
                {
                    checker = c ;
                }

                std::vector<tip>  TipDetect(const Graph & base )
                {
                    std::vector<tip> ret ;
                    std::set<NodeId> used;
                    for( auto & x : base.nodes )
                    {
                        if( used.find( x.id ) != used.end() )
                        {
                            continue ;
                        }
                        if( x.EdgeNum() == 1 )
                        {
                            tip tmp_tip ;
                            NodeId curr_id = x.id ;
                            tmp_tip.push_back(x.id);
                            const auto & node = base.GetNode(curr_id) ;
                            auto edge_id = *(node.edge_ids.begin());
                            const auto & edge = base.GetEdge(edge_id);
                            auto next = edge.OppoNode(curr_id) ;
                            while( checker(tmp_tip) )
                            {
                                const auto & node = base.GetNode(next) ;
                                if ( node.EdgeNum() == 0 )
                                {
                                    assert(0);
                                }
                                else if ( node.EdgeNum() == 1 )
                                {
                                    auto edge_id = *(node.edge_ids.begin());
                                    const auto & edge = base.GetEdge(edge_id);
                                    tmp_tip.push_back(next);
                                    next = edge.OppoNode(next) ;
                                }
                                else
                                { // tip end
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
                        return ret ;
                    }
                }

                void TipRemove( Graph & base , const tip & t)
                {
                    for(auto x : t )
                    {
                        base.RemoveNode(x);
                    }
                }

                //////////////////////////////////////////////////////
                //
                // Brief :
                //          cyclic remove tips from a graph.
                //
                //          loop :
                //              1. it is detect all tip .
                //              2. delete tips that shorter than tip_cut_num
                //
                /////////////////////////////////////////////////////
                TipRemoveResult DeepTipRemove( Graph & base )
                {
                    TipRemoveResult ret ;
                    ret.Init() ;
                    ret.base_node_num = base.nodes.size();
                    do
                    {
                        tip tips ;
                        tips = TipDetect( base );
                        if( tips.empty() )
                        {
                            break;
                        }
                        for( auto x : tips )
                        {
                            ret.tip_num ++ ;
                            ret.tip_node_num += x.size() ;
                            TipRemove(x);
                        }

                    }while( true );
                    ret.base_left_node_num = ret.base_node_num - ret.tip_node_num ;
                    return ret ;
                }
            };

    }
}


#endif //__ALGORITHM_GRAPH_TIP_REMOVE_H__
