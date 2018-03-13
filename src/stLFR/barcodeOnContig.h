#ifndef __STLFR_BARCODEONCONTIG_H__
#define __STLFR_BARCODEONCONTIG_H__

#include <set>
#include <map>
#include <vector>

namespace BGIQD {
namespace stLFR {

    struct MergedEdge
    {
        std::vector<unsigned int>  contigs;
    };

    struct SubGraphEdge
    {
        unsigned int id;
        std::set<unsigned int> subtos;
    };

    struct GraphEA_B
    {
        std::map<unsigned int , std::map<int,int> > barcode_on_contig;
    };

}// namespace stLFR
}// namespace BGIQD

#endif //__STLFR_BARCODEONCONTIG_H__
