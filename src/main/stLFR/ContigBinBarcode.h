#ifndef __STLFR_CBB_H__
#define __STLFR_CBB_H__

#include "utils/collection/collection.h"
#include "utils/incr_array/incr_array.h"
#include "utils/misc/contigIndex.h"
#include <map>
#include <set>

/**********************************************************
 *
 * @Brief : 
 *  Define some common data structures for store contig/bin/barcde
 *  and their relationships.
 *
 *  More than 50% of application in this scaffolder need this
 *  structures.
 *
 * *******************************************************/

namespace BGIQD{
    namespace stLFR {

        // The barcoding info for a contig
        // Store the which barcode aligned on this contig at where.
        struct ContigBarcodeInfo
        {
            unsigned int contig_id ;

            std::map<int , std::set<unsigned int > > barcodesOnPos;

            std::string ToString() const ;

            void Touch( int pos , unsigned int barcode)
            {
                barcodesOnPos[pos].insert(barcode);
            }

            void InitFromString(const std::string & line ) ;
        };

        // Store all contigs that have overlap with one barcode.
        // Use this index to avoid zero-Jaccard-similarity calculation.
        struct ContigOnBarcode
        {
            int barcode_id ;

            std::map<unsigned int , int > contig_data;

            void InitFromString(const std::string & line) ;

            std::string ToString() const ;

            void Touch(unsigned int contig , int num = 1)
            {
                if( contig_data.find(contig) == contig_data.end())
                    contig_data[contig] =num ;
                else
                    contig_data[contig]+=num ;
            }

        };

        // A set of barcodes. atomic item to calculate Jaccard similarty.
        typedef BGIQD::Collection::Collection<int> BarcodeCollection;

        // A bin with it's information : 
        //    + contig information
        //    + barcode information
        struct BarcodeOnBin
        {
            BGIQD::SOAP2::ContigId contigId ;

            int binId ;

            int start ;

            int end ;

            BarcodeCollection collections ;

            std::string ToString() const ;

            bool empty() const ;

            void InitFromString( const std::string & line ) ;
        };

        // The JS for 2 bin.
        // Only log the other bin's information because a bin always know itself.
        struct BinSimularity
        {
            int binIndex ;
            unsigned int contigId ;
            int binId ;
            float simularity;
        };

        // All JS for one bin
        struct BinRelation
        {
            int binIndex ; // bin index in all bins

            unsigned int contigId ;

            int binId ; // bin id in contig

            std::map<unsigned int , BinSimularity> sims;

            int start ;

            int end ;

            std::string ToString() const ;

            bool empty() const ;

            void InitFromString( const std::string & line ) ;
        };

        // The JS for 2 contig.
        // Only log the other contig because an contig always know itself.
        struct ContigSimularity
        {
            unsigned int contigId ;
            float simularity ;
        };

        // All JS for one contig
        struct ContigRelation
        {
            unsigned int contigId ;
            std::map<unsigned int , ContigSimularity> sims ;

            std::string ToString() const ;

            bool empty() const  { return sims.empty() ; }

            void InitFromString( const std::string & line ) ;
        };

        //
        // 
        // Below define some data matrix and their save and load functions :
        //
        //

        typedef  BGIQD::INCRARRAY::IncrArray<BarcodeOnBin>      BarcodeOnBinArray;
        void LoadBarcodeOnBinArray( const std::string & file , BarcodeOnBinArray & data );
        void PrintBarcodeOnBinArray( const std::string & file , const BarcodeOnBinArray & data);

        typedef  BGIQD::INCRARRAY::IncrArray<BinRelation>       BinRelationArray;
        void LoadBinRelationArray(const std::string & file , BinRelationArray & data);
        void PrintBinRelationArray(const std::string & file ,const BinRelationArray & data);

        typedef  BGIQD::INCRARRAY::IncrArray<ContigRelation>    ContigRelationArray;
        void LoadContigRelationArray(const std::string & file ,ContigRelationArray & data);
        void PrintContigRelationArray(const std::string & file ,const ContigRelationArray & data);

    }
}

#endif //__STLFR_CBB_H__
