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
            std::map<int , BinSimularity> sims;
        };

        typedef  BGIQD::INCRARRAY::IncrArray<BGIQD::stLFR::BarcodeOnBin>  BarcodeOnBinArray;

        void LoadBarcodeOnBinArray( const std::string & file , BarcodeOnBinArray & data );

        void PrintBarcodeOnBinArray( const std::string & file , const BarcodeOnBinArray & data);
    }
}

#endif //__STLFR_CBB_H__
