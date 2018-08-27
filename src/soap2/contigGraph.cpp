#include "soap2/contigGraph.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

#include <sstream>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace BGIQD{
    namespace SOAP2{

        // -------------------------- struct Edge ---------------------------------
        int Edge::ArcNum() const 
        {
            int ret = 0 ;
            Arc * cur_arc =  arc ;
            while( cur_arc != NULL )
            {
                ret ++ ; 
                cur_arc = cur_arc->next;
            }
            return ret ;
        }
/*
        void Edge::CheckLinear( Edge & a , Edge & b_a )
        {
            if( a.arc && (a.arc->next == 0) && b_a.arc && (b_a.arc->next == 0) )
            {
                a.flag |= 0x8 ;
                b_a.flag |= 0x8;
            }
        }

        void Edge::CheckTip( Edge &a , Edge &b_a )
        {
            if( a.arc && b_a.arc == 0 )
            {
                a.flag |= 0x10 ;
                b_a.flag |= 0x20;
            }
            else if ( a.arc == 0 && b_a.arc )
            {
                a.flag |= 0x20 ;
                b_a.flag |= 0x10;
            }
        }

        void Edge::CheckRepeat(Edge &a , Edge &b_a )
        {
            if( ( a.arc &&a.arc->next ) || ( b_a.arc && b_a.arc->next ) )
            {
                a.flag |= 0x2 ;
                b_a.flag |= 0x2;
            }
        }
*/
        void Edge::DepthSearch(Edge * array 
                ,std::list<Edge> & stack
                ,std::map<unsigned int , Edge> & history
                ,std::map<unsigned int , std::vector<std::list<Edge>> > & paths
                ,std::map<unsigned int , std::vector<std::list<Edge>> > & mids 
                ,int total_length
                ,const std::map<unsigned int , float> & neibs
                ,int max_depth )
        {
            if ( total_length >= max_depth|| length == 0 )
                return ;
            stack.push_front(*this);
            auto & top = stack.front();
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
                        stack.push_front(next);
                        paths[next.id].push_back(stack);
                        stack.pop_front();
                        step = true;
                        continue;
                    }
                    if( neibs.find( next.bal_id) != neibs.end() )
                    {
                        stack.push_front(next);
                        paths[next.bal_id].push_back(stack);
                        stack.pop_front();
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
                next.DepthSearch(array,
                        stack
                        ,history
                        ,paths
                        ,mids
                        ,total_length + next.length
                        ,neibs
                        ,max_depth);
            }
            if(step)
                top.JumpStep();
            history[top.id]=top;
            stack.pop_front();
        }

        // -------------------------- struct KeyEdge ------------------------------

        KeyConn & KeyEdge::GetValidTo() 
        {
            for( auto & i : to)
            {
                if( i.second.IsValid() )
                    return i.second;
            }
            assert(0);
            static KeyConn i ;
            return i;
        }
        KeyConn & KeyEdge::GetValidFrom() 
        {
            for( auto & i : from )
            {
                if( i.second.IsValid() )
                    return i.second;
            }
            assert(0);
            static KeyConn i ;
            return i;
        }
        void KeyConn::InitFromString( const std::string & str )
        {
            auto items = BGIQD::STRING::split(str,":");
            assert( items.size() == 4 ) ;
            to = std::stoul( items[0] );
            length = std::stoul( items[1] );
            sim = std::stof( items[2] );
            flag = 0 ;
            if( items[3] == "+" )
                SetPostive() ;
        }

        std::tuple<bool,bool,bool> KeyEdge::Relationship(unsigned int id) const 
        {
            auto itr1 = from.find(id);
            if( itr1 != from.end() )
            {
                return std::make_tuple(true , false , itr1->second.IsPositive());
            }

            auto itr2 = to.find(id);
            if( itr2 != to.end() )
            {
                return std::make_tuple(true , true , itr2->second.IsPositive() );
            }
            return std::make_tuple(false ,false ,false);
        }

        void KeyEdge::Init(int i , unsigned int b)
        {
            id = i ;
            edge_id = b ;
            bal_id = b + 1;
        }

        void KeyEdge::InitFrom( const std::string & str )
        {
            KeyConn conn ;
            conn.InitFromString( str );
            if ( conn.IsPositive() )
            {
                from[conn.to] = conn ;
            }
            else
            {
                from[conn.to + 1 ] = conn ;
            }
        }

        void KeyEdge::InitTo( const std::string & str )
        {
            KeyConn conn ;
            conn.InitFromString( str );
            if ( conn.IsPositive() )
            {
                to[conn.to] = conn ;
            }
            else
            {
                to[conn.to + 1 ] = conn ;
            }
        }

        std::tuple<bool,bool,bool> KeyEdge::Relationship_nojump(unsigned int id,bool to_order) const 
        {
            if( to_order )
            {
                auto itr2 = to.find(id);
                if(itr2 != to.end() && (! itr2->second.IsJumpConn()) )
                {
                    return std::make_tuple(true , true , itr2->second.IsPositive() );
                }
            }
            else
            {
                auto itr1 = from.find(id);
                if(itr1 != from.end() && (! itr1->second.IsJumpConn()))
                {
                    return std::make_tuple(true , false , itr1->second.IsPositive());
                }
            }
            return std::make_tuple(false ,false ,false);
        }

        void KeyEdge::CheckCircle()
        {
            for( const auto & f : from )
            {
                if( ! f.second.IsValid() )
                    continue ;
                if( to.find(f.first) != to.end() && to.at(f.first).IsValid() )
                {
                    flag |= 0x40;
                    return ;
                }
            }
            flag &= 0xffffffbf ;
        }

        void KeyEdge::SetType() 
        {
            flag = 0 ;
            from_size = 0 ; to_size = 0;
            total_size = 0 ;
            jump_conn = 0 ;
            for( const auto & i : from)
            {
                if( i.second.IsJumpConn() )
                {
                    jump_conn ++ ;
                    continue;
                }
                if( i.second.IsBiNotSupport() )
                {
                    jump_conn ++ ;
                    continue;
                }
                from_size ++ ;
            }
            for( const auto & i : to)
            {
                if( i.second.IsJumpConn() )
                {
                    jump_conn ++ ;
                    continue;
                }
                if( i.second.IsBiNotSupport() )
                {
                    jump_conn ++ ;
                    continue;
                }
                to_size ++ ;
            }
            total_size = from_size + to_size ;
            if( from_size == 1 && to_size == 1 )
                flag |= 0x1;
            else if( from_size == 0 && to_size == 0 )
                flag |= 0x20;
            else if ( from_size > 0 && to_size == 0 )
                flag |= 0x2;
            else if ( from_size == 0 && to_size > 0 )
                flag |= 0x10 ;
            else
                flag |= 0x4;
            CheckCircle();
        }

        // -------------------------- struct GraphEA------------------------------

        void GraphEA::LoadEdge( const std::string & file, int K)
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( "open prefix.update_edge for read failed !!! " );
            std::string line;
            std::getline(*in,line);
            sscanf(line.c_str(),"EDGEs %u",&contigTotalNum);
            edge_array = static_cast<Edge*>( calloc(sizeof(Edge) , contigTotalNum + 1 ));
            int length;
            int bal;
            int cov;
            unsigned int index = 1 ;
            while(!std::getline(*in,line).eof())
            {
                auto & from =  edge_array[index].from ;
                auto & to =  edge_array[index].to;
#if K127mer
                sscanf(line.c_str(),">length %d,%d,%d,%lx %lx %lx %lx,%lx %lx %lx %lx",
                        &length,&bal,&cov
                        ,&(from[0]),&(from[1]) ,&(from[2]),&(from[3])
                        ,&(to[0]),&(to[1]) ,&(to[2]),&(to[3]) );
#else
                sscanf(line.c_str(),">length %d,%d,%d,%lx %lx,%lx %lx",
                        &length,&bal,&cov
                        ,&(from[0]),&(from[1])
                        ,&(to[0]),&(to[1]));
#endif
                edge_array[index].id = index ;
                edge_array[index].bal_id = index+bal;
                edge_array[index].cov = cov / 10;
                if( length > K )
                    edge_array[index].length = length ;
                else
                    edge_array[index].length = 0 ;
                if( bal == 0 )
                    edge_array[index].SetPalindrome();
                if( bal == 1 )
                    edge_array[index].SetBase();

                if( edge_array[index].length == 0)
                {
                    edge_array[index].SetDelete();
                }
                index ++ ;
            }
            assert( index == contigTotalNum +1 );
            return ;
        }

        void GraphEA::LoadArc( const std::string & file)
        {
            std::string line;
            unsigned int contigId;
            unsigned int to;
            int cov;
            arcNum = 0;
            // Counting arcs
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( "open prefix.Arc for read failed !!! " );
            while(!std::getline(*in,line).eof())
            {
                std::istringstream ist(line);
                ist>>contigId;
                while(! ist.eof() )
                {
                    ist>>to>>cov;
                    arcNum ++ ;
                }
            }
            delete in ;
            arc_array =static_cast<Arc*>( calloc(sizeof(Arc),arcNum + 1));
            long long index = 1 ;
            auto in1 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( "open prefix.Arc for read failed !!! " );
            while(!std::getline(*in1,line).eof())
            {
                std::istringstream ist(line);
                ist>>contigId;
                while( ! ist.eof() )
                {
                    // make sure are list are sorted by to id .
                    ist>>to>>cov;
                    arc_array[index].to = to;
                    arc_array[index].cov = cov;

                    const auto & to_node = edge_array[to];
                    if( to_node.IsDelete() && to_node.length < 1 )
                    {
                        index++ ;
                        continue ;
                    }
                    auto & node = edge_array[contigId] ;

                    if (node.arc == NULL ||  to <=  node.arc->to )
                    {
                        arc_array[index].next = edge_array[contigId].arc;
                        edge_array[contigId].arc = &arc_array[index];
                    }
                    else
                    {
                        Arc ** prev = &node.arc ;
                        Arc * curr = node.arc ;
                        while( curr != NULL )
                        {
                            if( curr->to <= to ) 
                                break ;
                            prev = &(curr->next );
                            curr = curr->next ;
                        }
                        arc_array[index].next = curr ;
                        *prev = &arc_array[index];
                    }
                    index ++ ;
                }
            }
            delete in1;
            assert(index == arcNum +1 );
        }
    }//namespace SOAP2
}//namespace BGIQD
