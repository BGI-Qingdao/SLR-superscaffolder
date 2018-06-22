#ifndef __STLFR_CBB_H__
#define __STLFR_CBB_H__

#include "algorithm/collection/collection.h"
#include "algorithm/incr_array/incr_array.h"
#include "soap2/soap2.h"
#include <map>
#include <set>

//
// file stLFR/CBB.h
//
// struct of Contig &  Bin  & Barcode ...
//
//

namespace BGIQD{
    namespace stLFR {

        typedef BGIQD::Collection::Collection<int> BarcodeCollection;

        struct BarcodeOnBin
        {
            BGIQD::SOAP2::ContigId contigId ;

            int binId ;

            BarcodeCollection collections ;

            std::string ToString() const ;

            bool empty() const ;

            void InitFromString( const std::string & line ) ;
        };

        struct BinSimularity
        {
            int binIndex ;
            unsigned int contigId ;
            int binId ;
            float simularity;
        };

        struct BinRelation
        {
            int binIndex ;
            unsigned int contigId ;
            int binId ;
            std::map<unsigned int , BinSimularity> sims;

            std::string ToString() const ;

            bool empty() const ;

            void InitFromString( const std::string & line ) ;
        };

        struct ContigSimularity
        {
            unsigned int contigId ;
            float simularity ;
            //TODO : length , weight ...
        };
        struct ContigRelation
        {
            unsigned int contigId ;
            std::map<unsigned int , ContigSimularity> sims ;

            std::string ToString() const ;

            bool empty() const  { return sims.empty() ; }

            void InitFromString( const std::string & line ) ;
        };

        typedef  BGIQD::INCRARRAY::IncrArray<BarcodeOnBin>      BarcodeOnBinArray;
        typedef  BGIQD::INCRARRAY::IncrArray<BinRelation>       BinRelationArray;
        typedef  BGIQD::INCRARRAY::IncrArray<ContigRelation>    ContigRelationArray;

        void LoadBarcodeOnBinArray( const std::string & file , BarcodeOnBinArray & data );
        void PrintBarcodeOnBinArray( const std::string & file , const BarcodeOnBinArray & data);
        void LoadBinRelationArray(const std::string & file , BinRelationArray & data);
        void PrintBinRelationArray(const std::string & file ,const BinRelationArray & data);
        void LoadContigRelationArray(const std::string & file ,ContigRelationArray & data);
        void PrintContigRelationArray(const std::string & file ,const ContigRelationArray & data);

    }
}

#endif //__STLFR_CBB_H__
