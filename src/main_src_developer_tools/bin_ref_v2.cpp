#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
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
        if( pos < start ) return ;
        auto & bin = getByPos(pos);
        bin.data.IncreaseElement(barcode);
    }
    private:
    Bin & getByPos( int pos ) {
        return bins.at((pos-start )/ bin_size) ;
    }
};

struct AppConfig {

    std::vector<BinPos> NZoo;
    BinRef from ;
    BinRef to ;

    void Init( int ref_len , int bin_size , int dist ) {
        from.Init(0,bin_size,ref_len);
        to.Init(bin_size+dist,bin_size,ref_len);
    }

    void Loadread2contig(const std::string & file ) {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(file);
        if( NULL == in ) FATAL("failed to open read2contig file to read !");
        std::string line ;
        while( ! std::getline(*in, line).eof() ) {
            BGIQD::EASY_SAM::EasySam tmp ;
            tmp.InitFromString(line);
            from.assign(tmp.pos_1bp,tmp.barcode);
            to.assign(tmp.pos_1bp,tmp.barcode);
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
        from.MarkN(NZoo);
        to.MarkN(NZoo);
    }

    void ClacAndPrint() {
        std::cout<<"A_bt\tB_bt\tAUB_bt\tANB_bt\tA_bn\tB_bn\tAUB_bn\tANB_bn\n";
        int i ; int jump = 0 ;
        for( i = 0 ; i < from.size() ; i ++ ) {
            if(!from.has(i) || ! to.has(i)) { jump++ ; continue ; }
            const auto & bin_from = from.get(i);
            const auto & bin_to   =    to.get(i);
            auto AUB = BinData::Union(bin_from.data,bin_to.data);
            auto ANB = BinData::Intersection(bin_from.data,bin_to.data);
            std::cout<<bin_from.data.keysize()
                <<'\t'<<bin_to.data.keysize()
                <<'\t'<<AUB.keysize()
                <<'\t'<<ANB.keysize()
                <<'\t'<<bin_from.data.size()
                <<'\t'<<bin_to.data.size()
                <<'\t'<<AUB.size()
                <<'\t'<<ANB.size()<<'\n';
        }
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
        DEFINE_ARG_REQUIRED( int , bin_size , " the bin size ");
        DEFINE_ARG_REQUIRED( int , dist, " the dist between bin ");
        DEFINE_ARG_REQUIRED( int , ref_len, " the length of ref seq ");
    END_PARSE_ARGS;

    config.Init(ref_len.to_int() , bin_size.to_int() , dist.to_int()) ;
    config.Loadread2contig(read2contig.to_string());
    config.MaskNZoo(bed.to_string());
    config.ClacAndPrint();
    return 0;
}
