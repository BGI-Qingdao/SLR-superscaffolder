#include "stLFR/TagId.h"
#include <iostream>

namespace BGIQD{
    namespace stLFR{

        long long TagId::AddTag(const std::string & tag)
        {
            auto itr = m_tag2num.find(tag);
            if( itr != m_tag2num.end() )
                return itr->second;
            curr ++ ;
            m_tag2num[tag] = curr ;
            curr ++ ;
            return curr - 1;
        }

        long long TagId::Id(const std::string & tag ) const 
        {
            auto itr = m_tag2num.find(tag);
            if( itr != m_tag2num.end() )
                return itr->second;
            return -1;
        }

        void TagId::AssignTag( const std::string &tag,long long number )
        {
            m_tag2num[tag] = number;
        }

        void TagId::Print( std::ostream & out ) const 
        {
            for( const auto & i : m_tag2num )
            {
                out<<i.first<<'\t'<<i.second<<'\n';
            }
        }

        /*********************************************************************/
    }
}
