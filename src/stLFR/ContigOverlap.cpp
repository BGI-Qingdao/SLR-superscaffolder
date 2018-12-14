#include "stLFR/ContigOverlap.h"
#include <cassert>

namespace BGIQD {
    namespace stLFR {

        void OverlapInfo::InitFromPAF( const BGIQD::PAF::PAF_Item &item )
        {
            unsigned int contig1 = std::stoi(item.query_name);
            int c1len = item.query_len;
            int c1_start = item.query_start;
            int c1_end = item.query_end;
            char c1pos = item.query_char;
            unsigned int contig2 = std::stoi( item.target_name) ;
            int c2len  = item.target_len;
            int c2_start = item.target_start;
            int c2_end = item.target_end ;
            int len_m1 = item.len_query_match;
            //int len_m2 = item.len_target_match ;

            if( c2_start < c2len - c2_end )
            { // c1 --> c2
                if(  c1pos == '+' )
                {
                    if( c1_start < c1len - c1_end )
                    {
                        //     -------- <---c1+
                        //           -------- c2 -->
                        type = OverlapType::Unknow;
                    }
                    else
                    {
                        // +c1-->  ------
                        //           -------- c2 ->
                        type = OverlapType::C1_2_C2 ;
                    }
                }
                else
                {
                    assert( c1pos == '-' );
                    if( c1_start < c1len - c1_end )
                    {
                        //     -------- <---c1-
                        //           -------- c2 -->
                        type = OverlapType::C1p_2_C2 ;
                    }
                    else
                    {
                        // -c1--> ------
                        //           -------- c2 ->
                        type = OverlapType::Unknow;
                    }
                }
            }
            else 
            { // c2 --> c1
                // c1 --> c2p
                if(  c1pos == '+' )
                {
                    if( c1_start < c1len - c1_end )
                    {
                        //              +c1--> -------- 
                        //           ------------- c2 -->
                        type = OverlapType::C1p_2_C2p;
                    }
                    else
                    {
                        //                   -------- <--c1+
                        //           ------------- c2 -->
                        type = OverlapType::Unknow;
                    }
                }
                else
                {
                    assert( c1pos == '-' );
                    if( c1_start < c1len - c1_end )
                    {
                        //             -c1 --------
                        //           -------- c2 -->
                        type = OverlapType::Unknow;
                    }
                    else
                    {
                        //               ------  c1 -
                        //           -------- c2 ->
                        type = OverlapType::C1_2_C2p;
                    }
                }
            }
            c1 = contig1 ;
            c2 = contig2 ;

            if( c1 > c2 && type != OverlapType::Unknow )

            {
                std::swap(c1,c2);
                type = static_cast<OverlapType>(5-static_cast<int>(type)) ;
            }
            overlap_len = len_m1 ;
        }

        void OverlapInfo::InitFromScaffItem( const BGIQD::stLFR::ContigDetail & prev 
                , const BGIQD::stLFR::ContigDetail & next)
        {
            c1 = prev.contig_id ;
            c2 = next.contig_id ;

            if ( prev.orientation && next.orientation )
            {
                type = OverlapType::C1_2_C2 ;
            }
            else if ( prev.orientation && ( ! next.orientation ) )
            {
                type = OverlapType::C1_2_C2p ;
            }
            else if ( ( ! prev.orientation ) && (  next.orientation ) )
            {
                type = OverlapType::C1p_2_C2 ;
            }
            else
            {
                type = OverlapType::C1p_2_C2p ;
            }
            if( c1 > c2 && type != OverlapType::Unknow )

            {
                std::swap(c1,c2);
                type = static_cast<OverlapType>(5-static_cast<int>(type)) ;
            }
        }
    }
}
