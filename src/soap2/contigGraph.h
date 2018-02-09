#ifndef __SOAP2_CONTIGGRAPH_H__
#define __SOAP2_CONTIGGRAPH_H__

#include <vector>
#include <map>
#include <stack>
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
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > paths
            ,std::map<unsigned int , std::vector<std::stack<Edge>> > mids 
            ,int total_length
            ,const std::map<unsigned int , float> & neibs);
};

struct Connection
{
    virtual bool IsConnected () { return false ; }
    virtual bool IsUniqueConnected () { return false ; }
    virtual int Level1Count() { return 0 ; }
    virtual int TotalCount() { return 0 ;}
    unsigned int to ;
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

struct KeyEdge
{
    unsigned int id ;
    unsigned int edge_id ;
    int flag ;
    std::map<unsigned int ,Connection * > connections;
};

}
}

#endif //__SOAP2_CONTIGGRAPH_H__
