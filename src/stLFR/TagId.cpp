#include "stLFR/TagId.h"
#include <iostream>

namespace BGIQD{
    namespace stLFR{

        int TagId::AddTag(const std::string & tag)
        {
            auto itr = m_tag2num.find(tag);
            if( itr != m_tag2num.end() )
                return itr->second;
            while(1)
            {
                if( m_num2tag.find(curr) == m_num2tag.end() )
                    break;
                curr ++ ;
            }
            m_num2tag[curr]=tag;
            m_tag2num[tag] = curr ;
            curr ++ ;
            return curr - 1;
        }

        int TagId::Id(const std::string & tag ) const 
        {
            auto itr = m_tag2num.find(tag);
            if( itr != m_tag2num.end() )
                return itr->second;
            return -1;
        }

        bool TagId::AssignTag( const std::string &tag, int number )
        {
            auto itr = m_tag2num.find(tag);
            if( itr != m_tag2num.end() )
                return itr->second == number ;
            auto itr1 = m_num2tag.find(number);
            if( itr1 != m_num2tag.end() )
                return false ;
            m_num2tag[number]=tag;
            m_tag2num[tag] = number;
            return true ;
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
