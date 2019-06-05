#include "biocommon/agp/agp.h"

#include <sstream>

namespace BGIQD {
    namespace AGP {

        std::vector<std::string> split_by_tab(
                const std::string & line)
        {
            std::vector<std::string> ret ;
            std::string tmp="";
            for( char x : line )
            {
                if( x == '\t')
                {
                    if(tmp.size() > 0 )
                        ret.push_back(tmp);
                    tmp="";
                    continue ;
                }
                if( x == '\n' ) break ;
                tmp+=x;
            }
            if(tmp.size() > 0 )
                ret.push_back(tmp);

            return ret ;
        }

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

        void Scaff2AGPItem::InitName(const std::string & n)
        {
            scaff_name = n ;
            length = 1 ;
            part_number = 0 ;
        }

        void Scaff2AGPItem::AddSeq(
                const std::string & sname ,
                int sbegin ,
                int send , 
                char orientation
                )
        {
            long long add_length = send - sbegin +1 ;
            AGP_Item tmp ;
            tmp.object = scaff_name ;
            tmp.object_beg = length ;
            tmp.object_end = length + add_length -1;
            tmp.part_number = ++part_number ;
            tmp.component_type = AGP_Item::ComponentType::W ;
            tmp.lefta.component_id = sname ;
            tmp.lefta.orientation = orientation ;
            tmp.lefta.component_beg = sbegin ;
            tmp.lefta.component_end = send ;

            length += add_length ;
            data.emplace_back(std::move(tmp));
        }

        void Scaff2AGPItem::AddN(int n_size)
        {
            AGP_Item tmp ;
            tmp.object = scaff_name ;
            tmp.object_beg = length;
            tmp.object_end = length + n_size -1;
            tmp.part_number = ++part_number ;
            tmp.component_type = AGP_Item::ComponentType::N ;
            tmp.leftb.gap_length = n_size ;
            tmp.leftb.gap_type = "scaffold";
            tmp.leftb.linkage = true ;
            tmp.leftb.linkage_evidence = "map";

            length += n_size;
            data.emplace_back(std::move(tmp));
        }

        bool AGP_Item::IsAGPStringValid(const std::string &str) 
        {
            int tab_num = 0 ;
            for( const char c : str ) if ( c == '\t' ) tab_num ++ ;
            return tab_num == 8 ;
        }
        void AGP_Item::InitFromString(const std::string & str)
        {
//scaffold_1      1996387 1996396 172     N       10      scaffold        yes     map
//scaffold_1      1996397 2012814 173     W       170_264 1       16418   -
            char component_type_c ;
            std::istringstream ist(str);
            ist>>object>>object_beg>>object_end>>part_number>>component_type_c;
            component_type = static_cast<ComponentType>(component_type_c);
            if( component_type == ComponentType::N ||
                    component_type == ComponentType::U )
            {
                std::string linkage_str ;
                ist>>leftb.gap_length>>leftb.gap_type
                    >>linkage_str>>leftb.linkage_evidence;
                if( linkage_str == "yes")
                    leftb.linkage = true ;
                else 
                    leftb.linkage = false ;
            }
            else
            {
                ist>>lefta.component_id>>lefta.component_beg
                    >>lefta.component_end>>lefta.orientation ;
            }
        }
    }
}
