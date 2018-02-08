#ifndef __SOAP2_CONTIGGRAPH_H__
#define __SOAP2_CONTIGGRAPH_H__

#include <vector>
#include <map>
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

    bool IsDelete() const { return flag & 0x1 ; }
    bool IsRepeat() const { return flag & 0x2 ; }
    bool IsUnique() const { return flag & 0x4 ; }
    bool IsLinear() const { return flag & 0x8 ; }
    bool IsTipStart() const { return flag & 0x10 ; }
    bool IsTipEnd() const { return flag & 0x20 ; }

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
