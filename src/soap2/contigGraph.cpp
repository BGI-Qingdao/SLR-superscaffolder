#include "soap2/contigGraph.h"
#include <common/files/file_reader.h>
#include <sstream>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

namespace BGIQD{
    namespace SOAP2{


        // -------------------------- struct Edge ---------------------------------

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

        void Edge::DepthSearch(Edge * array 
                ,std::list<Edge> & stack
                ,std::map<unsigned int , Edge> & history
                ,std::map<unsigned int , std::vector<std::list<Edge>> > & paths
                ,std::map<unsigned int , std::vector<std::list<Edge>> > & mids 
                ,int total_length
                ,const std::map<unsigned int , float> & neibs)
        {
            if ( total_length >= 1000000 || length == 0 )
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
                        ,neibs);
            }
            if(step)
                top.JumpStep();
            history[top.id]=top;
            stack.pop_front();
        }

        // -------------------------- struct KeyEdge ------------------------------
        std::tuple<bool,bool,bool> KeyEdge::Relationship(unsigned int id)
        {
            auto itr1 = from.find(id);
            if(itr1 != from.end() )
            {
                return std::make_tuple(true , false , itr1->second.IsPositive());
            }

            auto itr2 = to.find(id);
            if(itr2 != to.end() )
            {
                return std::make_tuple(true , true , itr2->second.IsPositive() );
            }
            return std::make_tuple(false ,false ,false);
        }


        void KeyEdge::CheckCircle()
        {
            for( const auto & f : from )
            {
                if( to.find(f.first) != to.end() )
                {
                    flag |= 0x40;
                    return ;
                }
            }
        }

        void KeyEdge::SetType() 
        {
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
                from_size ++ ;
            }
            for( const auto & i : to)
            {
                if( i.second.IsJumpConn() )
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
                flag |= 0x10 ;
            else if ( from_size == 0 && to_size > 0 )
                flag |= 0x2;
            else
                flag |= 0x4;
            CheckCircle();
        }

        // -------------------------- struct GraphEA------------------------------

        void GraphEA::LoadEdge( const std::string & file, int K)
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
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
                sscanf(line.c_str(),">length %d,%d,%d",&length,&bal,&cov);
                edge_array[index].id = index ;
                edge_array[index].bal_id = index+bal;
                edge_array[index].cov = cov;
                if( length > K )
                    edge_array[index].length = length;
                else
                    edge_array[index].length = 0 ;
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
            while(!std::getline(*in1,line).eof())
            {
                std::istringstream ist(line);
                ist>>contigId;
                while(! ist.eof() )
                {
                    ist>>to>>cov;
                    arc_array[index].to = to;
                    arc_array[index].cov = cov;
                    arc_array[index].next = edge_array[contigId].arc;
                    edge_array[contigId].arc = &arc_array[index];
                    index ++ ;
                }
            }
            delete in1;
            assert(index == arcNum +1 );
        }

    }//namespace SOAP2
}//namespace BGIQD
