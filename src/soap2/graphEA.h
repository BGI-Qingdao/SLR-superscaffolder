#ifndef __SOAP2_GRAPHEA_H__
#define __SOAP2_GRAPHEA_H__

#include "soap2/kmer.h"
#include "algorithm/incr_array/incr_array.h"
#include "common/flags/flags.h"

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
            std::string     K;
            std::string     linear;

            int ArcNum() const ;

            FLAGS_INT;

            ADD_A_FLAG(1,Delete);
            ADD_A_FLAG(2,Key);
            ADD_A_FLAG(3,Base);
            ADD_A_FLAG(4,Mark);
            ADD_A_FLAG(5,Parlindorme);
            ADD_A_FLAG(6,Unique);
            ADD_A_FLAG(7,Repeat);
            ADD_A_FLAG(8,TipStart);
            ADD_A_FLAG(9,UsedInSuper);
            ADD_A_FLAG(10,UsedInLinear);
            ADD_A_FLAG(11,JumpStep);
            ADD_A_FLAG(12,IsTipEnd);
            ADD_A_FLAG(13,Merge);
            ADD_A_FLAG(14,SetK);


            void Init(const std::string & line , int k) ;

            Edge ReverseCompelete() const ;

            void AddSeq ( const std::string & line, int k) ;

            std::string ToString() const ;

            std::string ToString(int new_id, const std::string &marker ) const ;

            bool IsSeqComplete(int k) const
            {
                return ( (int)K.size() == k && (int)linear.size() ==length );
            }

        };

        // contig as vertex.
        // arc between contig as directed edge.
        // The Edge Arc graph of SOAPdenovo contig.
        //
        struct  GraphEA
        {
            typedef BGIQD::INCRARRAY::IncrArray<Edge>  EdgeArray;
            EdgeArray edge_array;
            typedef BGIQD::INCRARRAY::IncrArray<Arc> ArcArray;
            ArcArray arc_array;
            void LoadEdge( const std::string & file, int K);
            void LoadArc( const std::string & file);
            void LoadSeq(const std::string & file);
            void BuildReverseCompleteSeq();
        }; // struct GraphEA


    } // SOAP2
} // BGIQD

#endif //__SOAP2_GRAPHGEA_H__
