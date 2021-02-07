#ifndef __BIOCOMMON_AGP_H__
#define __BIOCOMMON_AGP_H__
/******************************************
 *
 * codes for print info in AGP format: "AGP Specification v2.0"
 * https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
 *
 * ****************************************/

#include <string>
#include <vector>

namespace BGIQD{

    // Define the classes for save scaffolding information by AGP format.
    namespace AGP {

        // Define one row of AGP as AGP_Item
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
            // If this is a contig row.
            struct PartA
            {
                std::string component_id;   // column 6a
                long long component_beg ;   // column 7a
                long long component_end ;   // column 8a
                std::string orientation ;   // column 9a
            };
            // If this is a gap row.
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
        };

        // Define rows of AGP_Items as AGPFile, and implement Print function.
        struct AGPFile
        {
            std::vector<AGP_Item> data;
            // implement.
            void Print( std::ostream & ost ) const ;
        };
    }
}

#endif //__BIOCOMMON_AGP_H__
