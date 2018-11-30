#ifndef __STLFT_TAGID_H__
#define __STLFT_TAGID_H__

#include <string>
#include <map>

namespace BGIQD {
    namespace stLFR {

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
                long AddTag(const std::string & tag );
                // Return the number of a tag , -1 for not found .
                long  Id( const std::string & tag ) const ;

                // Assign a number for a tag .
                // Return false if assign fail.
                //      like tag already exsist and number not match.
                //          or number alread exsist and tag not match.
                bool AssignTag( const std::string & tag , long number );

                // Print barcodeList into file
                //  format is :
                //      "tag\tid\n"
                void Print(std::ostream & ost ) const;
            private:
                std::map<std::string , long > m_tag2num;
                std::map<long , std::string > m_num2tag;
                long curr ;
        };

    }
}

#endif
