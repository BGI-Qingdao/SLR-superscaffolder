#include "biocommon/fasta/fasta.h"

#include <string>
#include <sstream>
namespace BGIQD {

    namespace FASTA {

        void SOAP2ContigHead::Init( const std::string & line ) 
        {
            //>21 length 64 cvg_0.0_tip_0
            sscanf(line.c_str() 
                    ,">%u length %d cvg_%f_tip_%d"
                    ,&contigId
                    ,&len
                    ,&cov 
                    ,&is_tip);

        }

        std::string SOAP2ContigHead::Head() const {
                std::ostringstream ost;
                ost<<'<'<<contigId
                    <<" length "<<len
                    <<" cvg_"<<cov
                    <<"_tip_"<<is_tip;
                return ost.str();
            }

    }
}
