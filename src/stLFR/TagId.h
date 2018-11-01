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
                int AddTag(const std::string & tag );
                // Return the number of a tag , -1 for not found .
                int  Id( const std::string & tag ) const ;

                // Assign a number for a tag .
                // Return false if assign fail.
                //      like tag already exsist and number not match.
                //          or number alread exsist and tag not match.
                bool AssignTag( const std::string & tag , int number );

                // Print barcodeList into file
                //  format is :
                //      "tag\tid\n"
                void Print(std::ostream & ost ) const;
            private:
                std::map<std::string , int > m_tag2num;
                std::map<int , std::string > m_num2tag;
                int curr ;
        };

    }
}

#endif
