#ifndef __STLFR_CONTIGOVERLAP_H__
#define __STLFR_CONTIGOVERLAP_H__

#include <string>
#include "stLFR/ScaffInfo.h"
#include "biocommon/paf/PAF.h"

namespace BGIQD {
    namespace stLFR {

        //
        // Any overlap can display like below :
        //
        //                  c1-start-index
        //                      |
        // contig c1 (c1/c1p)   V
        //      -------------------              contig c2 (c2/c2p)
        //      ^               ---------------------------
        //      |                 ^                       ^
        //  1-index-c1            |                       |
        //                        c2-end-index            1-index-c2
        //

        struct OverlapInfo
        {
            unsigned int c1 ;
            unsigned int c2 ;

            enum OverlapType
            {
                Unknow      = 0 ,
                C1_2_C2     = 1 ,
                C1_2_C2p    = 2 ,
                C1p_2_C2    = 3 ,
                C1p_2_C2p   = 4 ,
            };

            static inline OverlapType GetCR( OverlapType type )
            {
                if(  type == C1_2_C2 )
                {
                    return C1p_2_C2p ;
                }
                if( type == C1p_2_C2p )
                {
                    return C1_2_C2 ;
                }
                return type ; 
                //return static_cast<OOType>( 5- static_cast<int>(type) );
            }

            OverlapType type ;
            int overlap_len ;
            int c1_start_pos ;
            int c2_end_pos ;

            std::string Key() const 
            {
                return std::to_string(c1) +"_"+ std::to_string(c2);
            }

            std::string KeyT() const 
            {
                return std::to_string(c1) +"_"+ std::to_string(c2)+std::to_string(static_cast<int>(type));
            }

            void InitFromPAF( const BGIQD::PAF::PAF_Item & item );

            void InitFromScaffItem( const BGIQD::stLFR::ContigDetail & prev 
                    , const BGIQD::stLFR::ContigDetail & next);
        };
    }
}

#endif
