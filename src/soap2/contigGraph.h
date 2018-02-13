#ifndef __SOAP2_CONTIGGRAPH_H__
#define __SOAP2_CONTIGGRAPH_H__

#include <vector>
#include <map>
#include <stack>
#include <set>
#include <functional>
namespace BGIQD {
namespace SOAP2 {

struct Arc
{
    unsigned int to;
    int cov;
    Arc * next ;
};

struct Edge
{
    unsigned int id ;
    unsigned int bal_id ;

    int length;
    int cov;
    int flag ; // bits marker
    Arc * arc;

    //Flag 
    bool IsDelete() const { return flag & 0x1 ; }
    bool IsRepeat() const { return flag & 0x2 ; }
    bool IsUnique() const { return flag & 0x4 ; }
    bool IsLinear() const { return flag & 0x8 ; }
    bool IsTipStart() const { return flag & 0x10 ; }
    bool IsTipEnd() const { return flag & 0x20 ; }
    bool IsKey() const { return flag & 0x40 ; }
    bool IsJumpStep() const { return flag & 0x80 ; }

    void SetKey() { flag |= 0x40 ; }
    void JumpStep() { flag |= 0x80 ;}

    static void CheckLinear( Edge & a , Edge & b_a )
    {
        if( a.arc && (a.arc->next == 0) && b_a.arc && (b_a.arc->next == 0) )
        {
            a.flag |= 0x8 ;
            b_a.flag |= 0x8;
        }
    }

    static void CheckTip( Edge &a , Edge &b_a )
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

    static void CheckRepeat(Edge &a , Edge &b_a )
    {
        if( ( a.arc &&a.arc->next ) || ( b_a.arc && b_a.arc->next ) )
        {
            a.flag |= 0x2 ;
            b_a.flag |= 0x2;
        }
    }

    void DepthSearch(Edge * array 
            ,std::stack<Edge> & stack
            ,std::map<unsigned int , Edge> & history
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > &paths
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > &mids 
            ,int total_length
            ,const std::map<unsigned int , float> & neibs);
};

struct Connection
{
    virtual bool IsConnected () { return false ; }
    virtual bool IsUniqueConnected () { return false ; }
    virtual int Level1Count() { return 0 ; }
    virtual int TotalCount() { return 0 ;}
};

struct LinearConnection  : public Connection 
{
    std::vector<int> line;
};

struct BubbleConnection : public Connection
{
    std::vector<std::vector<int>> bubbles;
};

struct BubblesConnection : public Connection
{
    std::vector<std::vector<std::vector<int>>> all;
};

struct KeyConn
{
    unsigned int to;
    int length;
    int flag ;

    bool IsPositive() const { return flag & 0x2 ; }
    void SetPostive() { flag |= 0x2 ;}
    bool IsJumpConn() const { return flag & 0x1 ;}
    void SetJump() { flag |= 0x1 ; }
};

struct KeyEdge
{
    unsigned int id ;
    unsigned int edge_id ;
    int flag ;
    //std::map<unsigned int ,Connection * > connections;
    //  for each edge , it has positive and reverse order.
    //  let the small id be positice
    //  let the bigger id be reverse
    //  from map is the downstream of reverse order
    //  to map is the downstream of positive order
    std::map<unsigned int , KeyConn> from ;
    std::map<unsigned int , KeyConn> to;

    // connect ? up/down ? p/r
    std::tuple<bool,bool,bool> Relationship(unsigned int id)
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

    bool IsLinear() const { return flag & 0x1 ; }
    bool IsTipFrom() const { return flag & 0x2 ; }
    bool IsJunction() const { return flag & 0x4 ; }
    bool IsMarked() { return flag & 0x8 ;}
    bool IsTipTo() const { return flag & 0x10 ; }
    bool IsSingle() const { return flag & 0x20 ; }

    void Mark() { flag |= 0x8 ; }
    void SetType() 
    {
        int from_size = 0 , to_size = 0;
        for( const auto & i : from)
        {
            if( i.second.IsJumpConn() )
                    continue;
            from_size ++ ;
        }
        for( const auto & i : to)
        {
            if( i.second.IsJumpConn() )
                    continue;
            to_size ++ ;
        }

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
    }
};

}
}

#endif //__SOAP2_CONTIGGRAPH_H__
