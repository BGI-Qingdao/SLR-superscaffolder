#ifndef __STLFT_TAGID_H__
#define __STLFT_TAGID_H__

#include <string>
#include <map>

/**********************************************************
 *
 *  A string <---> number mapping system.
 *
 *  Use number to replace string to save MEMORY
 *
 * ********************************************************/
namespace BGIQD {
    namespace MISC {

        //
        //  Tag <---> Number 
        //      1       :       1
        //
        class TagId
        {
            public:
                TagId() : curr(1) {}
                TagId(const TagId & );
                TagId & operator =(const TagId & );
            public:
                // Add a tag and return it's number.
                long long AddTag(const std::string & tag );
                // Return the number of a tag , -1 for not found .
                long long  Id( const std::string & tag ) const ;

                // Assign a number for a tag .
                void AssignTag( const std::string & tag ,long long number );

                // Print barcodeList into file
                //  format is :
                //      "tag\tid\n"
                void Print(std::ostream & ost ) const;
            private:
                std::map<std::string ,long  long > m_tag2num;
                //std::map<long , std::string > m_num2tag;
                long curr ;
        };

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

                bool HasId( long id ) ;

                std::map<long , std::string > data;
        };
    }
}

#endif
