#include "biocommon/seq/seq.h"
#include <sstream>
namespace BGIQD {
    namespace SEQ {

            std::string seq::Seq(int weight ) const
            {
                if( weight < 1 )
                    return atcgs ;
                else
                {
                    std::ostringstream ost;
                    int i = 1 ;
                    for( char c : atcgs)
                    {
                        ost<<c;
                        if( i % weight == 0 || i == (int)atcgs.size() )
                            ost<<'\n';
                        i++ ;
                    }
                    return ost.str();
                }
            }
    }
}
