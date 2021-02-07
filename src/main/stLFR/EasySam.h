#ifndef __STLFR_EASY_SAM_H__
#define __STLFR_EASY_SAM_H__

#include <string>
#include <sstream>
#include <cassert>

/**********************************************************
 *
 * @Brief :
 *       define some structure to store alignments.
 *       used in 
 *          +Sam2ReadOnContig
 *          +SplitInfo
 *          +PEGraph
 *
 *********************************************************/

namespace BGIQD {

    namespace EASY_SAM
    {
        // info about a single aligment.
        struct EasySam
        {
            long long read_id ;
            unsigned int contig_name ;
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

        // Item for PEInfo.
        struct PE_Single
        {
            long long  read1;
            bool match_reverse1;
            unsigned int contig1;
            int pos_1bp1;

            std::string ToString() const 
            {
                std::ostringstream ost;
                ost<<read1<<'\t';
                ost<<contig1<<'\t';
                ost<<pos_1bp1<<'\t';
                ost<<match_reverse1;
                return ost.str();
            }

            void InitFromString(const std::string & line)
            {
                std::istringstream ist(line);
                ist>>read1
                    >>contig1
                    >>pos_1bp1
                    >>match_reverse1;
            }
        };

        // info for a proper mapped read pair. 
        struct PEInfo
        {
            long long read1;
            long long read2;
            bool match_reverse1;
            bool match_reverse2;
            unsigned int contig1;
            unsigned int contig2;
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

            PE_Single PInfo() const 
            {
                PE_Single tmp;
                tmp.read1 = read1 ;
                tmp.match_reverse1 = match_reverse1 ;
                tmp.pos_1bp1 = pos_1bp1 ;
                tmp.contig1 = contig1 ;
                return tmp ;
            };

            PE_Single EInfo() const 
            {
                PE_Single tmp;
                tmp.read1 = read2 ;
                tmp.match_reverse1 = match_reverse2 ;
                tmp.pos_1bp1 = pos_1bp2 ;
                tmp.contig1 = contig2 ;
                return tmp ;
            }
        };

        // A statistics for all read-pair
        struct PE_Baisc
        {
            int expect ;
            int sd ;
            int total_pair_in_same;
            int total_pair_in_diff;
        };
    }
}

#endif //__STLFR_EASY_SAM_H__
