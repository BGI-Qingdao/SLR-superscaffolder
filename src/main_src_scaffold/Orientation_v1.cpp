#include "common/args/argsparser.h"
#include "common/error/Error.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/freq/freq.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/string/stringtools.h"
#include "common/stl/mapHelper.h"
#include "common/stl/setHelper.h"

#include "soap2/fileName.h"
#include "soap2/contigIndex.h"
#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <map>
#include <vector>
#include <array>
#include <algorithm>

struct ContigBarcodeSet{
    unsigned int id ;
    std::set<int> left ;
    std::set<int> right ;
    std::set<int> all;
    bool valid() const {
        return  ! ( left.empty() || right.empty() || all.empty() ) ;
    }
};

std::map<unsigned int , ContigBarcodeSet> cbs ;

const ContigBarcodeSet & GetInfo( unsigned int id ) {
    static ContigBarcodeSet inv;
    if( cbs.find( id ) == cbs.end() ) {
        assert(0);
        return inv;
    }
    return cbs.at(id);
}

struct SupportInfo{
    int direction ;     // -1 ; 0 ; 1     ==> suppert left ; none ; suppert right
    int weight ;        //      0 ; 1     ==> none         ; support ;
    bool is_terminal;   //    true; false ==> terminal     ; not terminal
    void  Init(){
        weight = 0 ;
        direction = 0 ;
        is_terminal = false ;
    }
};

SupportInfo GetSupportInfo( 
        unsigned int neib 
        , unsigned int center
        , int min_shared 
        , float min_ration) {
    const auto & neib_info   = GetInfo(neib) ;
    const auto & center_info = GetInfo(center);
    SupportInfo ret ; ret.Init() ;
    if( ! neib_info.valid() || ! center_info.valid() ) return ret ;
    //   A   BC  ; use contig A to determine the direction of BC
    auto AB  = BGIQD::STL::set_common(neib_info.all , center_info.left);
    auto AC  = BGIQD::STL::set_common(neib_info.all , center_info.right);
    auto ABC = BGIQD::STL::set_common(AB , AC );
    if( (int)ABC.size() >= min_shared ) {
        auto AB1 = BGIQD::STL::set_diff_in_s1(center_info.left, neib_info.all);
        auto B2 = BGIQD::STL::set_diff_in_s1(center_info.left, center_info.right);
        float B2_ration = float(B2.size()) / float(center_info.all.size() ) ;
        auto AC1 = BGIQD::STL::set_diff_in_s1(center_info.right, neib_info.all);
        auto C2 = BGIQD::STL::set_diff_in_s1(center_info.right, center_info.left);
        float C2_ration = float(C2.size()) / float(center_info.all.size() ) ;
        if( B2_ration > C2_ration && B2_ration/C2_ration >= min_ration ) {
            ret.weight = 1 ;
            ret.direction = -1 ;
        }
        if( C2_ration >B2_ration && C2_ration/B2_ration >= min_ration ) {
            ret.weight = 1 ;
            ret.direction = 1 ;
        }
    }else {
        if( AB.size() > AC.size() ){
            ret.weight = 1 ;
            ret.direction = -1 ;
        }
        if( AB.size() < AC.size() ){
            ret.weight = 1 ;
            ret.direction = 1 ;
        }
        if( AB.empty() && AC.empty() )
            ret.is_terminal = true ;
    }
    return ret ;
}

struct AppConfig
{

    BGIQD::SOAP2::FileNames fName ;
    BGIQD::LOG::logger loger;

    struct  OrientationResult {
        unsigned int contig_id ;
        int direction ; // -1 ; 0 ; 1   ==>     -   ;   none    ;   +
        int weight_1  ; 
        int weight_n1 ; 
        int left_proven ; 
        int right_proven ;
        void Init(unsigned int id ) {
            contig_id = id ;
            direction = 0 ;
            weight_1 = 0 ;
            weight_n1 = 0 ;
            left_proven = 0 ;
            right_proven = 0 ;
        }
        void AddLeftProven(int d) {
            left_proven ++ ; 
            if( d == 1 ) weight_n1 ++ ;// left neib think right bin is closer ==> - orientation
            if( d == -1 ) weight_1 ++ ;// left neib think left  bin is closer ==> + orientation
        }
        void AddRightProven(int d) {
            right_proven ++ ; 
            if( d == 1 ) weight_1  ++;// right neib think right bin is closer ==> + orientation
            if( d == -1 ) weight_n1++;// right neib think left  bin is closer ==> - orientation
        }
        void GenVoteReuslt() { 
            if( weight_1 > weight_n1 ) direction = 1 ; // + orientation
            if( weight_1 < weight_n1 ) direction = -1 ;// - orientation 
            if( weight_1 == weight_n1) direction = 0 ;
        }
        std::string ToString() const {
            std::ostringstream ost ;
            ost<<contig_id<<'\t'<<direction<<'\t'<<weight_1<<'\t';
            ost<<weight_n1<<'\t'<<left_proven<<'\t'<<right_proven;
            return ost.str() ;
        }
    };
    std::map<unsigned int , OrientationResult> oresults ;
    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("Trunk2Scaff ",BGIQD::LOG::loglevel::INFO, loger);
    }
    std::map<unsigned int , std::vector<unsigned int> > scaffs;
    void LoadTrunk(){
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        unsigned int id = 0;
        std::string line ;
        while(! std::getline(*in,line).eof() )
        {
            if( line[0] == '-' )
            {
                id ++ ;
                continue ;
            }
            unsigned int now = std::stoul(line);
            scaffs[id].push_back(now);
            oresults[now].Init(now);
        }
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load Trunk done "<<BGIQD::LOG::lend() ;
    }

    BGIQD::SOAP2::ContigIndexMap contigIndexs ;
    void LoadContigIndexs() {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fName.ContigIndex());
        if(in == NULL)
            FATAL(" failed to open xxx.contigIndex for read!!! ");
        contigIndexs.LoadContigIndexs(*in);
        delete in ;
        contigIndexs.BuildReverseCompleteContigs();
    }

    void LoadBarcodeOnContig() {
        loger<<BGIQD::LOG::lstart() << "load barcodeOnContig start ..."<<BGIQD::LOG::lend() ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.BarcodeOnContig());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnContig to read !");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            int mid = contigIndexs.GetContigIndex(tmp.contig_id).length / 2 ;
            for( const auto & pair : tmp.barcodesOnPos) {
                for( int barcode : pair.second) {
                    int pos = pair.first ;
                    cbs[tmp.contig_id].all.insert(barcode);
                    if( pos < mid ) cbs[tmp.contig_id].left.insert(barcode);
                    else cbs[tmp.contig_id].right.insert(barcode);
                }
            }
        }
        delete in;
        loger<<BGIQD::LOG::lstart() << "load barcodeOnContig done"<<BGIQD::LOG::lend() ;
    }
    int min_shared ;
    float min_ration ;
    void VoteByNiebs( const std::vector<unsigned int > & scaff , int center_index ) {
        // check left neibs
        unsigned int center_id =  scaff.at(center_index) ;
        for( int i = center_index -1 ; i >= 0 ; i -- ) {
            auto check_ret = GetSupportInfo(scaff.at(i) , center_id, min_shared , min_ration );
            if( check_ret.weight == 1 ) 
                oresults[center_id].AddLeftProven(check_ret.direction) ;
            if( check_ret.is_terminal ) 
                break ;
        }
        // check right neibs
        for( int i = center_index -1 ; i < (int)scaff.size() ; i ++ ) {
            auto check_ret = GetSupportInfo(scaff.at(i) , center_id, min_shared , min_ration );
            if( check_ret.weight == 1 ) 
                oresults[center_id].AddRightProven(check_ret.direction) ;
            if( check_ret.is_terminal ) 
                break ;
        }
    };
    void GenOrientation() {
        for( const auto & pair : scaffs ) {
            for( int i = 0 ; i < (int)pair.second.size() ; i++ ) {
               VoteByNiebs( pair.second ,i ); 
            }
        }
        for( auto & pair : oresults ) {
            pair.second.GenVoteReuslt();
        }
    }
    void PrintOrientation() {
        auto out3 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.orientation());
        if( out3 == NULL )
            FATAL(" failed to open xxx.orientation for write !!! ");
        for( auto & i: oresults)
        {
            (*out3)<<i.second.ToString()<<'\n';
        }
        delete out3;
    }
} config;

int main(int argc, char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix of read names.\n\
                                                    In \n\
                                                        xxx.mintree_trunk_linear; xxx.barcodeOnContig ;\n\
                                                    Out\n\
                                                        xxx.gap_oo ;");
    END_PARSE_ARGS;

    config.Init( prefix.to_string());
    config.LoadContigIndexs();
    config.LoadTrunk();
    config.LoadBarcodeOnContig();
    config.GenOrientation();
    config.PrintOrientation();
    return 0;
}
