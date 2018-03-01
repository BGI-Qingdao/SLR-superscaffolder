#include "soap2/contigGraph.h"

namespace BGIQD{
namespace SOAP2{

    void Edge::DepthSearch(Edge * array 
            ,std::stack<Edge> & stack
            ,std::map<unsigned int , Edge> & history
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > & paths
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > & mids 
            ,int total_length
            ,const std::map<unsigned int , float> & neibs)
    {
        if( total_length >= 1000000 || length == 0 )
            return ;
        stack.push(*this);
        auto & top = stack.top();
        bool step = false;
        while(top.arc)
        {
            Edge & next = array[top.arc->to];
            Edge & next_bal = array[next.bal_id];
            top.arc = top.arc->next;
            if( next.IsKey() || next_bal.IsKey() )
            {
                if(neibs.find( next.id ) != neibs.end() ) 
                {
                    stack.push(next);
                    paths[next.id].push_back(stack);
                    stack.pop();
                    step = true;
                    continue;
                }
                if( neibs.find( next.bal_id) != neibs.end() )
                {
                    stack.push(next);
                    paths[next.bal_id].push_back(stack);
                    stack.pop();
                    step = true;
                    continue;
                }
                continue;
            }

            auto itr1 = history.find(next.bal_id);
            if( itr1 != history.end() )
            {
                if(itr1->second.IsJumpStep() )
                {
                    step = true;
                    mids[next.bal_id].push_back(stack);
                }
                continue;
            }
            auto itr = history.find(next.id);
            if( itr != history.end() )
            {
                if(itr->second.IsJumpStep() )
                {
                    step = true;
                    mids[next.id].push_back(stack);
                }
                continue;
            }
            next.DepthSearch(array,stack,history,paths,mids,total_length + next.length,neibs);
        }
        if(step)
            top.JumpStep();
        history[top.id]=top;
        stack.pop();
    }
}
}
