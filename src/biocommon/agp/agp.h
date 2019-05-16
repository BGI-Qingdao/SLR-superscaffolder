#ifndef __BIOCOMMON_AGP_H__
#define __BIOCOMMON_AGP_H__
/******************************************
 *
 * Tools code for "AGP Specification v2.0"
 * https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
 *
 * ****************************************/

#include <string>
#include <vector>

namespace BGIQD{
    namespace AGP {

        struct AGP_Item
        {
            std::string object;     //column 1
            long long object_beg ;  //column 2 
            long long object_end ;  //column 3 
            int  part_number;       //column 4
            enum ComponentType
            {
                A = 'A' ,
                D = 'D' ,
                F = 'F' ,
                G = 'G' ,
                O = 'O' ,
                P = 'P' ,
                W = 'W' ,
                N = 'N' ,
                U = 'U' ,
                Unknow = 0 
            }
            component_type ;        //column 5
            struct PartA
            {
                std::string component_id;   // column 6a
                long long component_beg ;   // column 7a
                long long component_end ;   // column 8a
                std::string orientation ;   // column 9a
            };
            struct PartB
            {
                int gap_length ;            // column 6b
                std::string gap_type;       // column 7b
                bool  linkage;              // column 8b
                std::string linkage_evidence; // column 9b
            };

            PartA lefta ;
            PartB leftb ;

            std::string ToString() const;

            void InitFromString() ;
        };

        struct AGPFile
        {
            std::vector<AGP_Item> data;

            void Load( std::istream & ist );
            void Print( std::ostream & ost ) const ;
        };

        struct Scaff2AGPItem
        {
            public:
                void InitName(const std::string & n);
                // sbegin && send is 1base index
                void AddSeq(const std::string & sname ,
                        int sbegin ,
                        int send , 
                        char orientation );
                void AddN(int n_size);

                const std::vector<AGP_Item> & Items() const { return data ; }
            private:
                long long length ;
                int part_number ;
                std::string scaff_name ;
                std::vector<AGP_Item> data;
        };
    }
}

#endif //__BIOCOMMON_AGP_H__
