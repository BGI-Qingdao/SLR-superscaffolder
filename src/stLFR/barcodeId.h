#ifndef __STLFR_BARCODEID_H__
#define __STLFR_BARCODEID_H__

#include <map>
#include <string>

namespace BGIQD{
namespace stLFR{

//
//  BarcodeTag <---> Number 
//      1       :       1
//  Not multi-thread safe !
//
class BarcodeId
{
    public:
        static BarcodeId & Singleton() { return the_one ; }
    private:
        BarcodeId() : curr(1) {}
        BarcodeId(const BarcodeId & );
        BarcodeId & operator =(const BarcodeId & );
        static BarcodeId the_one;

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

class BarcodeIdHelper
{
    public:
        // Is barcodeList preLoaded or not ?
        static bool preload ;
        // Return the id of tag 
        //  if preload == true
        //      return correct id or -1 
        //  else
        //      return correct id or assign a new id and return it .
        static int  Id( const std::string & tag );

        // Load barcodeList from file 
        //  content format :
        //      string\tint\n
        //  means:
        //      barcode\tid\n
        static void Load( const std::string & file);

        static void Print( const std::string & file);
};

} // namespace stLFR
} // namespace BGIQD
#endif //__STLFR_BARCODEID_H__
