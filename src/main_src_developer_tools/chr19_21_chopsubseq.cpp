#include "biocommon/fasta/fasta.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "biocommon/seq/seq.h"
#include "common/args/argsparser.h"
#include <cstdlib>
#include <ctime>
#include <sstream>
#include <random>
#include <chrono>
#include <algorithm>

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

FragmentInfo subfrag(  const BGIQD::SEQ::seq & base , int begin , int len ){
    FragmentInfo ret ;
    ret.from = begin  ;
    ret.len  = len ;
    ret.orient = true ;
    ret.seq.AddPartSeq(base.atcgs.substr(begin,len));
    return ret ;
}

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

int main(int argc , char ** argv ) {
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(int , bin_size, "bin size");
    END_PARSE_ARGS
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
    // gen correct seq in bin_size
    std::vector<FragmentInfo> correct_bin_19;
    std::vector<FragmentInfo> correct_bin_21;
    for( int i = 0 ; i < chr19.Len() ; i+=bin_size.to_int() ) {
        auto ret = subfrag(chr19 ,i,bin_size.to_int());
        ret.base = "chr19";
        correct_bin_19.push_back(ret);
    }
    for( int i = 0 ; i < chr21.Len() ; i+=bin_size.to_int() ) {
        auto ret = subfrag(chr21 ,i,bin_size.to_int());
        ret.base = "chr21";
        correct_bin_21.push_back(ret);
    }
    // gen merge mis-assmbly by bin_size/2 + bin_size/2
    std::vector<FragmentInfo> sub_bin_19;
    std::vector<FragmentInfo> sub_bin_21;
    for( int i = 0 ; i < chr19.Len() ; i+=bin_size.to_int()/2 ) {
        auto ret = subfrag(chr19 ,i,bin_size.to_int()/2);
        ret.base = "chr19";
        sub_bin_19.push_back(ret);
    }
    for( int i = 0 ; i < chr21.Len() ; i+=bin_size.to_int()/2 ) {
        auto ret = subfrag(chr21 ,i,bin_size.to_int()/2);
        ret.base = "chr21";
        sub_bin_21.push_back(ret);
    }
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(sub_bin_19.begin() , sub_bin_19.end() ,std::default_random_engine(seed));
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(sub_bin_21.begin() , sub_bin_21.end() ,std::default_random_engine(seed));
    std::vector<FragmentMergeInfo> mis_bin_chr19;
    std::vector<FragmentMergeInfo> mis_bin_chr21;
    for( int i = 0 ; i < (int)sub_bin_19.size()-1; i += 2 ){
        FragmentMergeInfo tmp ;
        tmp.left = sub_bin_19[i];
        tmp.right = sub_bin_19[i+1];
        mis_bin_chr19.push_back(tmp);
    }
    for( int i = 0 ; i < (int)sub_bin_21.size()-1; i += 2 ){
        FragmentMergeInfo tmp ;
        tmp.left = sub_bin_21[i];
        tmp.right = sub_bin_21[i+1];
        mis_bin_chr21.push_back(tmp);
    }
    // print all
    print("chr19_correct.fa",correct_bin_19);
    print("chr21_correct.fa",correct_bin_21);
    print("chr19_merge.fa",mis_bin_chr19);
    print("chr21_merge.fa",mis_bin_chr21);
    return 0 ;
}
