#ifndef __STLFR_TRUNKGAP_H__
#define __STLFR_TRUNKGAP_H__

#include <iostream>
#include <vector>

namespace BGIQD {
    namespace stLFR {

        template<class T>
            struct TrunkGap
            {
                typedef T  Data;
                unsigned int prev ;
                unsigned int next ;
                Data data;
            };

        template<class T>
            void Load_MST_Trunk_Linear( std::istream & ist , std::vector<TrunkGap<T> > & data )
            {
                typedef TrunkGap<T> DataItem;
                std::string line ;
                unsigned int prev = -1 ;
                while(! std::getline(ist,line).eof() )
                {
                    if( line[0] == '-' )
                    {
                        prev = -1 ;
                        continue ;
                    }
                    unsigned int now = std::stoul(line);
                    if( prev != (unsigned int )-1 )
                    {
                        DataItem info ;
                        info.prev = prev ;
                        info.next = now ;
                        data.push_back(info);
                    }
                    prev = now ;
                }
            }
    }
}

#endif //__STLFR_TRUNKGAP_H__
