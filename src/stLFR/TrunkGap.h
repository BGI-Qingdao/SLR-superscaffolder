#ifndef __STLFR_TRUNKGAP_H__
#define __STLFR_TRUNKGAP_H__

#include <iostream>
#include <sstream>
#include <vector>
#include <map>

namespace BGIQD {
    namespace stLFR {

        struct GapFill
        {
            unsigned int prev ;
            unsigned int next ;
            unsigned int true_prev;
            unsigned int true_next ;
            std::vector<unsigned int> extra;

            void InitFromString( const std::string & str)
            {
                std::istringstream ist(str); ;
                ist>>prev>>next>>true_prev>>true_next;
                unsigned int f ;
                while(! ist.eof() )
                {
                    ist>>f ;
                    extra.push_back(f);
                }
            }

            std::string ToString() const 
            {
                std::ostringstream ost ;
                ost<<prev<<'\t'<<next<<'\t'<<true_prev<<'\t'<<true_next;
                for(unsigned int i : extra)
                {
                    ost<<'\t'<<i;
                }
                return ost.str() ;
            }
        };
        template<class T>
            struct TrunkGap
            {
                typedef T  Data;
                unsigned int prev ;
                unsigned int next ;
                Data data;
            };

        template<class T>
            void Load_MST_Trunk_Linear( std::istream & ist , std::map<int ,std::vector<TrunkGap<T> > >& data )
            {
                typedef TrunkGap<T> DataItem;
                std::string line ;
                unsigned int prev = -1 ;
                int id = 0;
                while(! std::getline(ist,line).eof() )
                {
                    if( line[0] == '-' )
                    {
                        id ++ ;
                        prev = -1 ;
                        continue ;
                    }
                    unsigned int now = std::stoul(line);
                    if( prev != (unsigned int )-1 )
                    {
                        DataItem info ;
                        info.prev = prev ;
                        info.next = now ;
                        data[id].push_back(info);
                    }
                    prev = now ;
                }
            }
    }
}

#endif //__STLFR_TRUNKGAP_H__
