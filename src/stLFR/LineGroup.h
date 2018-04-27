#ifndef __STLFR_LINEGROUP_H__
#define __STLFR_LINEGROUP_H__
#include <vector>
#include <string>

namespace BGIQD {
    namespace stLFR {

        struct ContigRoad
        {
            bool headin;
            bool tailin;
            int linear_length;
            std::vector<unsigned int> group;

            std::string  to_string() const ;

            void init(const std::string & buff);

            bool needMerge() const { return linear_length > 1 ; };

            std::pair<unsigned int, unsigned int>  getLinearStep(int index);

            std::vector<unsigned int > contig_path;


            void AddGroup( const std::vector<unsigned int > & a );

            enum FillStatus
            {
                None = 0 ,
                Conflict = 1 ,
                PartSucc = 2 ,
                Complete = 3
            } status;

            int fill_num;
            std::vector<int> circle_runs;
        };

        struct ContigRoads
        {
            void LoadRoads( const std::string & file );
            std::vector<ContigRoad>  roads;
        };


        typedef std::vector<unsigned int> ContigRoadFill ;
        struct ContigRoadFills
        {
            void LoadContigRoadFills( const std::string & file );
            std::vector<ContigRoadFill>  fills;
            std::vector<bool> has_circle ; 
        };
    }
}

#endif //__STLFR_LINEGROUP_H__
