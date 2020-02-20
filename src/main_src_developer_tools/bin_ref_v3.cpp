#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

#include "algorithm/linear_fitting/Minimum_multiplication.h"
#include "algorithm/collection/collection.h"
#include "algorithm/statistics/common.h"
#include "algorithm/interval/Interval.h"

#include "stLFR/EasySam.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <tuple>

typedef BGIQD::Collection::Collection<int> BinData;
typedef BGIQD::INTERVAL::Interval<int> BinPos ;

struct Bin{
    BinPos pos ;
    BinData data ;
    bool   n_marker ;
};

struct BinRef {

    bool has(int bin_id) const {
        if( bins.find(bin_id) == bins.end() )
            return false ;
        return ! bins.at(bin_id).n_marker ;
    }

    const Bin & get( int bin_id ) const {
        assert( has(bin_id) );
        return bins.at(bin_id);
    }

    std::map<int , Bin>  bins;

    int start ;
    int bin_size ;
    int size() const { return bins.size() ; }

    void Init( int start_ , int bin_size_ , int total_len ) {
        start = start_ ;
        bin_size = bin_size_ ;
        int index = 0 ;
        for( int i = start ; i < total_len  ; i+= bin_size ) {
            bins[index].pos.min = i;
            bins[index].pos.max = i+bin_size;
            bins[index].n_marker = false ;
            index ++ ;
        }
    }

    void MarkN( const std::vector<BinPos> & ns ){
        for( auto & pair : bins ) {
            auto & data = pair.second ;
            int n = 0 ;
            for( const auto & i : ns ) {
                n+= (i.Overlap(data.pos)).Len();
            }
            if( n *100 / data.pos.Len()  > 1 ) 
                data.n_marker = true ;
        }
    }
    void assign( int pos , int barcode ) {
        auto & bin = getByPos(pos);
        if( pos < start ) return ;
        bin.data.IncreaseElement(barcode);
    }
    private:
    Bin & getByPos( int pos ) {
        return bins.at((pos-start )/ bin_size) ;
    }
};

struct AppConfig {

    std::vector<BinPos> NZoo;
    BinRef bin1s;
    BinRef bin2s;
    BinRef bin3s;

    void Init( int ref_len , int bin_size , int dist1 , int dist2 ) {
        bin1s.Init(0,bin_size,ref_len);
        bin2s.Init(bin_size+dist1,bin_size,ref_len);
        bin3s.Init(bin_size+dist1+bin_size+dist2,bin_size,ref_len);
    }

    void Loadread2contig(const std::string & file ) {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(file);
        if( NULL == in ) FATAL("failed to open read2contig file to read !");
        std::string line ;
        while( ! std::getline(*in, line).eof() ) {
            BGIQD::EASY_SAM::EasySam tmp ;
            tmp.InitFromString(line);
            bin1s.assign(tmp.pos_1bp,tmp.barcode);
            bin2s.assign(tmp.pos_1bp,tmp.barcode);
            bin3s.assign(tmp.pos_1bp,tmp.barcode);
        }
        delete in ;
    }

    void MaskNZoo(const std::string & file ) {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(file);
        if( NULL == in ) FATAL("failed to open NZoo file to read !");
        std::string line ;
        while( ! std::getline(*in, line).eof() ) {
            std::string name ; 
            std::istringstream ist(line);
            BinPos tmp ; 
            ist>>name>>tmp.min>>tmp.max;
            NZoo.push_back(tmp);
        }
        delete in ;
        bin1s.MarkN(NZoo);
        bin2s.MarkN(NZoo);
        bin3s.MarkN(NZoo);
    }

    void ClacAndPrint(const std::string & prefix) {
        auto bn = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(prefix+".bn.txt");
        auto bt = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(prefix+".bt.txt");
        (*bn)<<"A_only\tB_only\tC_only\tABnC\tBCnA\tACnB\tABC\tF23\tAB\tAC\tBC\n";
        (*bt)<<"A_only\tB_only\tC_only\tABnC\tBCnA\tACnB\tABC\tF23\tAB\tAC\tBC\n";
        int i ; int jump = 0 ;
        for( i = 0 ; i < bin1s.size() ; i ++ ) {
            if(!bin1s.has(i)
            || !bin2s.has(i)
            || !bin3s.has(i) )
            { jump++ ; continue ; }
            const auto & A= bin1s.get(i).data;
            const auto & B= bin2s.get(i).data;
            const auto & C= bin3s.get(i).data;
            auto A_only   = A.Diff(B).Diff(C);
            auto B_only   = B.Diff(A).Diff(C);
            auto C_only   = C.Diff(A).Diff(B);
            auto AB       = BinData::Intersection(A,B);
            auto AC       = BinData::Intersection(A,C);
            auto BC       = BinData::Intersection(B,C);
            auto ABC      = BinData::Intersection(AB ,BC ) ;
            auto ABnC     = AB.Diff(C);
            auto ACnB     = AC.Diff(B);
            auto BCnA     = BC.Diff(A);
            float F23_bn  = B_only.size() == 0 ? 0 : ( float(std::min(A_only.size() ,C_only.size()) ) / float(B_only.size()) );
            float F23_bt  = B_only.keysize() == 0 ? 0 : ( float(std::min(A_only.keysize() ,C_only.keysize()) ) / float(B_only.keysize()) );
            (*bn)<<A_only.size()
                <<'\t'<<B_only.size()
                <<'\t'<<C_only.size()
                <<'\t'<<ABnC.size()
                <<'\t'<<BCnA.size()
                <<'\t'<<ACnB.size()
                <<'\t'<<ABC.size()
                <<'\t'<<F23_bn
                <<'\t'<<AB.size()
                <<'\t'<<AC.size()
                <<'\t'<<BC.size()<<"\n";
            (*bt)<<A_only.keysize()
                <<'\t'<<B_only.keysize()
                <<'\t'<<C_only.keysize()
                <<'\t'<<ABnC.keysize()
                <<'\t'<<BCnA.keysize()
                <<'\t'<<ACnB.keysize()
                <<'\t'<<ABC.keysize()
                <<'\t'<<F23_bt
                <<'\t'<<AB.keysize()
                <<'\t'<<AC.keysize()
                <<'\t'<<BC.keysize()<<"\n";
        }
        delete bn ;
        delete bt ;
        std::cerr<<" jump "<<jump<<" within "<<i<<" bin pars \n";
    }
} config ;

/******************************************************************************
 *
 * Logic.
 *
 *****************************************************************************/

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , read2contig , "the read2contig file");
        DEFINE_ARG_REQUIRED(std::string , bed , "the ref bed of N zoo");
        DEFINE_ARG_REQUIRED(std::string , prefix , "the output prefix");
        DEFINE_ARG_REQUIRED( int , bin_size , " the bin size ");
        DEFINE_ARG_REQUIRED( int , dist1, " the dist between bin1 & bin2 ");
        DEFINE_ARG_REQUIRED( int , dist2, " the dist between bin2 & bin3 ");
        DEFINE_ARG_REQUIRED( int , ref_len, " the length of ref seq ");
    END_PARSE_ARGS;

    config.Init(ref_len.to_int() , bin_size.to_int() , dist1.to_int() , dist2.to_int());
    std::string name_pref = prefix.to_string() 
                          + ".bin_"+std::to_string(bin_size.to_int())
                          + ".dist1_"+std::to_string(dist1.to_int())
                          + ".dist2_"+std::to_string(dist2.to_int());
    config.Loadread2contig(read2contig.to_string());
    config.MaskNZoo(bed.to_string());
    config.ClacAndPrint(name_pref);
    return 0;
}
