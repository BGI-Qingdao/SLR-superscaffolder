#ifndef __SOAP2_CONTIGGRAPH_H__
#define __SOAP2_CONTIGGRAPH_H__

#include <vector>
#include <map>
#include <stack>
#include <list>
#include <set>
#include <functional>
#include "soap2/kmer.h"
namespace BGIQD {
    namespace SOAP2 {


        // --------------------Graph Edge Arc --------------------
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
            Kmer from ;
            Kmer to ;

            //Flag
            void SetDelete()        { flag |= 0x1 ; }
            void SetUnique()        { flag |= 0x4 ; }
            void SetLinear()        { flag |= 0x8 ; }
            void SetKey()           { flag |= 0x40 ; }
            void JumpStep()         { flag |= 0x80 ; }
            void SetPalindrome()    { flag |= 0x100 ; }
            void SetBase()          { flag |= 0x200 ; }
            void SetMarked()        { flag |= 0x400 ; }

            bool IsDelete() const   { return flag & 0x1 ; }
            bool IsRepeat() const   { return flag & 0x2 ; }
            bool IsUnique() const   { return flag & 0x4 ; }
            bool IsLinear() const   { return flag & 0x8 ; }
            bool IsTipStart() const { return flag & 0x10 ; }
            bool IsTipEnd() const   { return flag & 0x20 ; }
            bool IsKey() const      { return flag & 0x40 ; }
            bool IsJumpStep() const { return flag & 0x80 ; }
            bool IsPalindrome()const{ return flag & 0x100 ; }
            bool IsBase() const     { return flag & 0x200 ; }
            bool IsMarked() const   { return flag & 0x400 ; }
            int ArcNum() const ;
            bool LinearDetect() const ;
            /*
            static void CheckLinear( Edge & a , Edge & b_a );

            static void CheckTip( Edge &a , Edge &b_a );

            static void CheckRepeat(Edge &a , Edge &b_a );
            */
            // Depth seach with neibs
            void DepthSearch(Edge * array 
                    ,std::list<Edge> & stack
                    ,std::map<unsigned int , Edge> & history
                    ,std::map<unsigned int , std::vector<std::list<Edge>> > &paths
                    ,std::map<unsigned int , std::vector<std::list<Edge>> > &mids 
                    ,int total_length
                    ,const std::map<unsigned int , float> & neibs
                    ,int max_length = 1000000);

        };

        // contig as vertex.
        // arc between contig as directed edge.
        // The Edge Arc graph of SOAPdenovo contig.
        //
        struct  GraphEA
        {
            Edge * edge_array;
            Arc * arc_array;
            unsigned int contigTotalNum;
            long long arcNum;
            void LoadEdge( const std::string & file, int K);
            void LoadArc( const std::string & file);
        }; // struct GraphEA

        // --------------------Key Edge&Cluster --------------------
        struct KeyConn
        {
            unsigned int to;
            int length;
            int flag ;

            bool IsPositive() const { return flag & 0x2 ; }
            void SetPostive() { flag |= 0x2 ;}

            bool IsJumpConn() const { return flag & 0x1 ;}
            void SetJump() { flag |= 0x1 ; }

            void SetBiNotSuppert() { flag |= 0x4 ; }
            bool IsBiNotSupport() const { return flag & 0x4 ;}


            bool IsValid() const { return  ! ( IsBiNotSupport() ||  IsJumpConn() ); }
        };

        struct KeyEdge
        {
            unsigned int id ;
            unsigned int edge_id ;
            unsigned int bal_id;
            int flag ;
            //std::map<unsigned int ,Connection * > connections;
            //  for each edge , it has positive and reverse order.
            //  let the small id ( edge_id of this struct ) be positice
            //  let the bigger id ( bal_id of edge_id ) be reverse
            //  from map is the downstream of reverse order
            //  to map is the downstream of positive order
            std::map<unsigned int , KeyConn> from ;
            std::map<unsigned int , KeyConn> to;

            // connect ? up/down ? p/r
            std::tuple<bool,bool,bool> Relationship(unsigned int id) const ;
            std::tuple<bool,bool,bool> Relationship_nojump(unsigned int id, bool to_order) const ;

            bool IsLinear() const { return flag & 0x1 ; }
            bool IsTipFrom() const { return flag & 0x2 ; }
            bool IsJunction() const { return flag & 0x4 ; }
            bool IsMarked() { return flag & 0x8 ;}
            bool IsTipTo() const { return flag & 0x10 ; }
            bool IsSingle() const { return flag & 0x20 ; }
            bool IsCircle() const { return flag & 0x40 ; }

            void Mark() { flag |= 0x8 ; }

            void CheckCircle();

            void SetType() ;

            int from_size ;
            int to_size ;
            int total_size;
            int jump_conn;
        }; // struct KeyEdge
    }//namespace SOAP2
}//namespace BGIQD

#endif //__SOAP2_CONTIGGRAPH_H__
