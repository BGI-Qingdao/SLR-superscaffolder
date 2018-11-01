#ifndef __STLFR_EASY_SAM_H__
#define __STLFR_EASY_SAM_H__

#include <string>
#include <sstream>
#include <cassert>

namespace BGIQD {

    namespace EASY_SAM
    {
        struct EasySam
        {
            long read_id ;
            int contig_name ;
            int pos_1bp;
            int barcode ;
            bool match_reverse;
            bool is_p ;
            bool pe_match;
            int insert_size ;

            std::string ToString() const 
            {
                std::ostringstream ost;
                ost<<read_id<<'\t';
                ost<<contig_name<<'\t';
                ost<<pos_1bp<<'\t';
                ost<<match_reverse<<'\t';
                ost<<barcode<<'\t';
                ost<<is_p<<'\t';
                ost<<pe_match<<'\t';
                ost<<insert_size;
                return ost.str();
            }

            void InitFromString(const std::string & line)
            {
                std::istringstream ist(line);
                ist>>read_id
                    >>contig_name
                    >>pos_1bp
                    >>match_reverse
                    >>barcode
                    >>is_p
                    >>pe_match
                    >>insert_size;
            }

        };

        struct PEInfo
        {
            long  read1;
            long  read2;
            bool match_reverse1;
            bool match_reverse2;
            int contig1;
            int contig2;
            int pos_1bp1;
            int pos_1bp2;

            std::string ToString() const 
            {
                std::ostringstream ost;
                ost<<read1<<'\t';
                ost<<contig1<<'\t';
                ost<<pos_1bp1<<'\t';
                ost<<match_reverse1<<'\t';
                ost<<read2<<'\t';
                ost<<contig2<<'\t';
                ost<<pos_1bp2<<'\t';
                ost<<match_reverse2;
                return ost.str();
            }

            void InitFromString(const std::string & line)
            {
                std::istringstream ist(line);
                ist>>read1
                    >>contig1
                    >>pos_1bp1
                    >>match_reverse1
                    >>read2
                    >>contig2
                    >>pos_1bp2
                    >>match_reverse2;
            }
        };

        struct PE_Baisc
        {
            int expect ;
            int sd ;
            int total_pair_in_same;
            int total_pair_in_diff;
        };

        struct EasySam_V1
        {
            long read_id ;

            int read_index;

            int contig_name ;

            int left_1bp;

            int barcode ;

            bool match_reverse;

            std::string ToString() const 
            {
                std::ostringstream ost;
                ost<<read_id<<'\t';
                ost<<read_index<<'\t';
                ost<<contig_name<<'\t';
                ost<<match_reverse<<'\t';
                ost<<left_1bp<<'\t';
                ost<<barcode;
                return ost.str();
            }

            void InitFromString(const std::string & line)
            {
                std::istringstream ist(line);
                ist>>read_id
                    >>read_index
                    >>contig_name
                    >>match_reverse
                    >>left_1bp
                    >>barcode;
            }
        };
    }
}

#endif //__STLFR_EASY_SAM_H__
