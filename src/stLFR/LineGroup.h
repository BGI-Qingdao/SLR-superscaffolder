#ifndef __STLFR_LINEGROUP_H__
#define __STLFR_LINEGROUP_H__

#include <vector>
#include <string>
#include "common/flags/flags.h"

namespace BGIQD {
    namespace stLFR {

        /**************** For fillcontigroad below *********************/
        struct ContigRoad
        {
            std::string  to_string() const ;

            void init(const std::string & buff);

            void AddGroup( const std::vector<unsigned int > & a );

            bool needMerge() const { return linear_length > 1 ; };

            std::pair<unsigned int, unsigned int>  getLinearStep(int index);

            // linear group
            bool headin;

            bool tailin;

            int linear_length;

            std::vector<unsigned int> group;

            // fill result
            std::vector<unsigned int > contig_path;

            int fill_num;

            std::vector<int> circle_runs;

            enum FillStatus
            {
                None = 0 ,
                Conflict = 1 ,
                PartSucc = 2 ,
                Complete = 3
            } status;

            ContigRoad Left(int fill)
            {
                ContigRoad ret ;
                ret.headin = ret.tailin = true ;
                if( fill >= linear_length - 1 )
                    return ret ;
                ret.linear_length = linear_length - fill - 1 ;
                ret.group.insert(ret.group.end(), group.begin() + fill , group.end() );
                return ret ;
            }
        };

        struct ContigRoads
        {
            void LoadRoads( const std::string & file );
            std::vector<ContigRoad>  roads;
        };

        /**************** For mergecontig below *********************/

        struct fill_flag
        {
            FLAGS_INT ;
            ADD_A_FLAG(0 , circle);
            ADD_A_FLAG(1 , fill_by_super);
        };

        typedef std::vector<unsigned int> ContigRoadFill ;

        struct ContigRoadFills
        {
            void LoadContigRoadFills( const std::string & file );

            std::vector<ContigRoadFill>  fills;

            std::vector<fill_flag> flags;
        };
    }
}

#endif //__STLFR_LINEGROUP_H__
