#include "biocommon/fasta/fasta.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "biocommon/seq/seq.h"
#include <cstdlib>
#include <ctime>
#include <sstream>

typedef BGIQD::FASTA::Fasta<BGIQD::FASTA::NormalHead> FA;
typedef BGIQD::SEQ::seq Seq;

struct FragmentInfo {
    std::string base ;
    int from ; // 0 base 
    int len ;
    bool orient;
    BGIQD::SEQ::seq seq;

    std::string ToString() const { 
        std::ostringstream ost ;
        ost<<'>'<<base<<"_begin:"<<from<<"_len:"<<len<<"_orient:"<<orient<<'\n';
        ost<<seq.Seq(60);
        return ost.str();
    }
    std::string Head() const {
        std::ostringstream ost ;
        ost<<base<<"_begin:"<<from<<"_len:"<<len<<"_orient:"<<orient;
        return ost.str();
    }
};

struct FragmentMergeInfo {
    FragmentInfo left ;
    FragmentInfo right;
    std::string ToString() const { 
        std::ostringstream ost ;
        ost<<">left_"<<left.Head()<<"_right_"<<right.Head()<<'\n';
        BGIQD::SEQ::seq tmp ;
        if( ! left.orient ) 
            tmp.AddPartSeq(left.seq.ReverseCompleteSeq());
        else 
            tmp.AddPartSeq(left.seq.atcgs);
        if( ! right.orient ) 
            tmp.AddPartSeq(right.seq.ReverseCompleteSeq());
        else 
            tmp.AddPartSeq(right.seq.atcgs);
        ost<<tmp.Seq(60);
        return ost.str();
    }
};

FragmentInfo fragment( const BGIQD::SEQ::seq & base , int len_min , int len_max ) {
    int begin = std::rand() % base.Len() ;
    int len_ran = 0 ;
    if ( len_max - len_min > 1 ) len_ran = std::rand() %( len_max - len_min);
    int len = len_min + len_ran ;
    FragmentInfo ret ;
    ret.from = begin  ;
    ret.len  = len ;
    ret.orient = true ;
    ret.seq.AddPartSeq(base.atcgs.substr(begin,len));
    if( ret.seq.NLen() * 100 / ret.seq.Len() > 10 ) return fragment(base ,len_min,len_max) ;
    return ret ;
}

BGIQD::SEQ::seq subseq(const BGIQD::SEQ::seq & base , int begin , int size ) {
    BGIQD::SEQ::seq ret ;
    ret.atcgs = base.atcgs.substr(begin,size);
    return ret;
}

template <class T> 
    void print ( const std::string & file ,const std::vector<T> & buffer ) {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
        for( const auto & x : buffer) {
            (*out)<<x.ToString();
        }
        delete out ;
    }

int main() {
    std::srand(std::time(nullptr));
    //file pos
    std::string chr21_file="01/chr21.fa";
    std::string chr19_file="01/chr19.fa";
    // chromesome ref
    Seq chr21 ;
    Seq chr19 ;
    // load ref 19 
    auto chr19_ist = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(chr19_file);
    std::vector<FA> buffer;
    BGIQD::FASTA::FastaReader<FA> reader ;
    reader.LoadAllFasta( *chr19_ist ,  buffer);
    assert(buffer.size() == 1 );
    chr19 = buffer[0].seq ;
    delete chr19_ist; buffer.clear();
    // load ref 21
    auto chr21_ist = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(chr21_file);
    std::vector<FA> buffer1;
    BGIQD::FASTA::FastaReader<FA> reader1 ;
    reader1.LoadAllFasta( *chr21_ist ,  buffer1);
    assert(buffer1.size() == 1 );
    chr21 = buffer1[0].seq ;
    delete chr21_ist; buffer1.clear();
    // generator random 10 50K framents earch 
    std::vector<FragmentInfo> type_b_19 ;
    std::vector<FragmentInfo> type_b_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr19 ,49000,50000);
        ret.base = "chr19";
        type_b_19.push_back(ret);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr21 ,49000,50000);
        ret.base = "chr21";
        type_b_21.push_back(ret);
    }
    // generate random 10  500~10000 fragments earch
    std::vector<FragmentInfo> type_c_19 ;
    std::vector<FragmentInfo> type_c_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr19 ,500,10000);
        ret.base = "chr19";
        type_c_19.push_back(ret);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr21 ,500,10000);
        ret.base = "chr21";
        type_c_21.push_back(ret);
    }
    // generate mis-assemble fragments
        // -- short short
    std::vector<FragmentMergeInfo> type_d_ss_19 ;
    std::vector<FragmentMergeInfo> type_d_ss_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr19 ,500,10000);
        auto ret2 = fragment(chr19 ,500,10000);
        ret.base = "chr19";
        ret2.base = "chr19";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_ss_19.push_back(tmp);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr21 ,500,10000);
        auto ret2 = fragment(chr21 ,500,10000);
        ret.base = "chr21";
        ret2.base = "chr21";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_ss_21.push_back(tmp);
    }
        // -- long -short
    std::vector<FragmentMergeInfo> type_d_sl_19 ;
    std::vector<FragmentMergeInfo> type_d_sl_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr19 ,500,10000);
        auto ret2 = fragment(chr19 ,30000,40000);
        ret.base = "chr19";
        ret2.base = "chr19";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_sl_19.push_back(tmp);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr21 ,500,10000);
        auto ret2 = fragment(chr21 ,30000,40000);
        ret.base = "chr21";
        ret2.base = "chr21";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_sl_21.push_back(tmp);
    }
        // -- long -long 
    std::vector<FragmentMergeInfo> type_d_ll_19 ;
    std::vector<FragmentMergeInfo> type_d_ll_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret2 = fragment(chr19 ,30000,40000);
        auto ret = fragment(chr19 ,30000,40000);
        ret.base = "chr19";
        ret2.base = "chr19";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_ll_19.push_back(tmp);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret2 = fragment(chr21 ,30000,40000);
        auto ret = fragment(chr21 ,30000,40000);
        ret.base = "chr21";
        ret2.base = "chr21";
        FragmentMergeInfo tmp ;
        tmp.left = ret ;
        tmp.right = ret2 ;
        type_d_ll_21.push_back(tmp);
    }
        // -- inversion 
    std::vector<FragmentMergeInfo> type_e_19 ;
    std::vector<FragmentMergeInfo> type_e_21 ;
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr19 ,30000,50000);
        ret.base = "chr19";
        int cut = 0 ;
        while( cut < 1005 || ret.seq.Len() - cut +1 < 1005 ) {
            cut = std::rand() % ret.seq.Len() ;
        }
        FragmentInfo left = ret ;
        left.seq = subseq( ret.seq , 0 , cut ) ;
	left.len = cut ;

        FragmentInfo right= ret ;
	right.from += cut ;
	right.len -= cut ;
        right.seq = subseq( ret.seq , cut , ret.seq.Len() - cut +1 );

        FragmentMergeInfo tmp ;
        tmp.left = left;
        tmp.left.orient = false ;
        tmp.right = right;
        type_e_19.push_back(tmp);
    }
    for( int i = 0 ; i < 10 ; i++ ) {
        auto ret = fragment(chr21 ,30000,50000);
        ret.base = "chr21";
        int cut = 0 ;
        while( cut < 1005 || ret.seq.Len() - cut +1 < 1005 ) {
            cut = std::rand() % ret.seq.Len() ;
        }

        FragmentInfo left = ret ;
        left.seq = subseq( ret.seq , 0 , cut ) ;
	left.len = cut ;

        FragmentInfo right= ret ;
	right.from += cut ;
	right.len -= cut ;
        right.seq = subseq( ret.seq , cut , ret.seq.Len() - cut +1 );

        FragmentMergeInfo tmp ;
        tmp.left = left;
        tmp.left.orient = false ;
        tmp.right = right;
        type_e_21.push_back(tmp);
    }
    // print all
    print("chr19_type_b.fa",type_b_19);
    print("chr19_type_c.fa",type_c_19);
    print("chr19_type_d_ss.fa",type_d_ss_19);
    print("chr19_type_d_sl.fa",type_d_sl_19);
    print("chr19_type_d_ll.fa",type_d_ll_19);
    print("chr19_type_e.fa",type_e_19);

    print("chr21_type_b.fa",type_b_21);
    print("chr21_type_c.fa",type_c_21);
    print("chr21_type_d_ss.fa",type_d_ss_21);
    print("chr21_type_d_sl.fa",type_d_sl_21);
    print("chr21_type_d_ll.fa",type_d_ll_21);
    print("chr21_type_e.fa",type_e_21);
    return 0 ;
}
