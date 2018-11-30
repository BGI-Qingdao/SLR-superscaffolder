#ifndef __STLFR_STRINGIDCACHE_H__
#define __STLFR_STRINGIDCACHE_H__

#include <map>
#include <string>
#include "stLFR/TagId.h"

namespace BGIQD{
    namespace stLFR{

        class StringIdCache
        {
            public:

                StringIdCache() : preload (false) {}

                // Is barcodeList preLoaded or not ?
                bool preload ;
                // Return the id of tag 
                //  if preload == true
                //      return correct id or -1 
                //  else
                //      return correct id or assign a new id and return it .
                long  Id( const std::string & tag );

                // Load barcodeList from file 
                //  content format :
                //      string\tint\n
                //  means:
                //      barcode\tid\n
                void Load( const std::string & file);

                void Print( const std::string & file);

                TagId data;
        };

        class IdStringCache
        {
            public:

                IdStringCache() : preload (false) {}

                // Is barcodeList preLoaded or not ?
                bool preload ;

                std::string  Id(long id);

                void LoadStringIdCache( const std::string & file);

                std::map<long , std::string > data;
        };
    } // namespace stLFR
} // namespace BGIQD
#endif //__STLFR_STRINGIDCACHE_H__
