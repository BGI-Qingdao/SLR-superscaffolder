#ifndef __STLFR_SCAFFINFO_H__
#define __STLFR_SCAFFINFO_H__

#include <vector>
#include <iostream>
#include <map>

/**********************************************************
 *
 * @Brief :
 *  Define a ScaffInfo format to store scaffolding infomations.
 *
 *  This format was replaced by AGP at final result but not all
 *  codes have been re-writed.
 *
 * The ScaffInfo format was used in :
 *   + Trunk2ScaffInfo
 *   + ScaffInfo2Seq
 *
 * *******************************************************/

namespace BGIQD {
    namespace stLFR {

        // Define details of a in-scaffold contig with 
        struct ContigDetail
        {
            unsigned int contig_id ;
            bool orientation ;
            int gap_size ;
            int contig_len ;
            int start_pos ;
            int scaff_index ;
            int scaff_id ;
            void InitFromString(const std::string & line) ;

            std::string ToString() const ;
        };

        // A scaffold : a collection of contigs
        struct ScaffInfo
        {
            int scaff_id ;
            std::vector<ContigDetail> a_scaff ;

            void PrintScaff( std::ostream & ost ) const;

            void FormatIndex();
            void FormatStartPos();
        };

        // Handle scaffolds : many scaffold and their relation with contigs
        struct ScaffInfoHelper
        {
            void PrintAllScaff( std::ostream & ost ) const ;

            void LoadAllScaff( std::istream & ist ) ;

            std::map<int , ScaffInfo> all_scaff ;

            struct ContigIndex
            {
                int scaffold_id ;
                int contig_index ;
            } ;

            std::map<unsigned int , ContigIndex>  contig_indexs ;

            ContigDetail & GetContig( unsigned int contig_id ) ;

            void FormatAllIndex()
            {
                for( auto & pair : all_scaff )
                    pair.second.FormatIndex();
            }
            void FormatAllStartPos()
            {
                for( auto & pair : all_scaff )
                    pair.second.FormatStartPos();
            }
        };
    }
}
#endif 
