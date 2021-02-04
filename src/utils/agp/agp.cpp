#include "utils/agp/agp.h"

#include <sstream>

namespace BGIQD {
    namespace AGP {

        std::string AGP_Item::ToString() const 
        {
            std::ostringstream ost ;
            ost<<object<<'\t'
               <<object_beg<<'\t'
               <<object_end<<'\t'
               <<part_number<<'\t'
               <<char(component_type)<<'\t';

            if( component_type ==ComponentType::N 
            || component_type == ComponentType::U )
            {
                const PartB &b = leftb ;
                ost<<b.gap_length<<'\t'
                   <<b.gap_type<<'\t'
                   <<(b.linkage ? "yes":"no")<<'\t'
                   <<b.linkage_evidence;
            }
            else
            {
                const PartA &a = lefta ;
                ost<<a.component_id<<'\t'
                 <<a.component_beg<<'\t'
                 <<a.component_end<<'\t'
                 <<a.orientation ;
            }
            return ost.str();
        }

        void AGPFile::Print( std::ostream & ost ) const
        {
             time_t now = time(0);
             char* dt = ctime(&now);
            // print comment first
            ost<<"##agp-version	2.0\n# ORGANISM: \n# TAX_ID: \n# ASSEMBLY NAME: \n";
            ost<<"# ASSEMBLY DATE: "<<dt;//<<'\n'; //ctime string have a \n
            ost<<"# GENOME CENTER: \n# DESCRIPTION: \n";
            // print item line by line
            for( const auto & item : data )
            {
                ost<<item.ToString()<<'\n';
            }
        }
    }
}
