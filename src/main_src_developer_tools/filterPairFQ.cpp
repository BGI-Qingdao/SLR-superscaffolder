#include "biocommon/fastq/fastq.h"
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include <string>
#include <map>
#include <set>

typedef BGIQD::FASTQ::stLFRHeader Header;
typedef BGIQD::FASTQ::Fastq<Header> Read;
typedef BGIQD::FASTQ::FastqReader<Read> Reader;

std::map<std::string , Read> reads1;
std::map<std::string , Read> reads2;
int main(int argc , char ** argv )
{
    
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,input_r1,"the input read1");
        DEFINE_ARG_REQUIRED(std::string,input_r2,"the input read2");
        DEFINE_ARG_REQUIRED(std::string,output_r1,"the output read1");
        DEFINE_ARG_REQUIRED(std::string,output_r2,"the output read2");
        DEFINE_ARG_REQUIRED(std::string,output_rs,"the output singlenton reads");
    END_PARSE_ARGS

    auto i_r1 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(input_r1.to_string());
    auto i_r2 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(input_r2.to_string());
    auto o_r1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_r1.to_string());
    auto o_r2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_r2.to_string());
    auto o_rs = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output_rs.to_string());
    if( i_r1 == NULL )
        FATAL("input_r1 failed to open !!");
    if( ! o_r1 || !o_r2 || ! o_rs )
        FATAL("open output file failed !! ");
    if( i_r2 == NULL )
        FATAL("input_r2 failed to open !!");

    std::vector<Read> buffer;
    // Load read1
    Reader reader1;
    reader1.LoadAllFastq(*i_r1,buffer);
    std::cerr<<" Load "<<buffer.size()<<" reads from read1. "<<std::endl;
    for( const auto & r : buffer )
    {
        reads1[r.head.readName] = r ;
    }
    buffer.clear();
    delete i_r1 ;
    // Load read2
    Reader reader2;
    reader2.LoadAllFastq(*i_r2,buffer);
    std::cerr<<" Load "<<buffer.size()<<" reads from read2. "<<std::endl;
    for( const auto & r : buffer )
    {
        reads2[r.head.readName] = r ;
    }
    buffer.clear();
    delete i_r2 ;
    // find common ;
    std::set<std::string> common ;
    for( const auto & pair : reads1 )
    {
        if( reads2.find( pair.first ) != reads2.end() )
        {
            common.insert(pair.first) ;
        }
    }
    
    std::cerr<<" Total "<<common.size()<<" read pair detected. "<<std::endl;
    auto print_stlfr_reads = []( const Read & r , std::ostream & ost ) {
            ost<<r.head.Head()<<'\n'<<r.seq.atcgs<<"\n+\n"<<std::string(r.seq.atcgs.size() , 'F')<<'\n';
    };
    std::cerr<<"Start print read1 ... "<<std::endl;
    // print read1 
    for( const auto & pair : reads1 )
        if( common.find(pair.first) != common.end() )
            print_stlfr_reads(pair.second , *o_r1);
    delete o_r1 ;
    std::cerr<<"Start print read2 ... "<<std::endl;
    // print read2 
    for( const auto & pair : reads2 )
        if( common.find(pair.first) != common.end() )
            print_stlfr_reads(pair.second , *o_r2);
    delete o_r2 ;
    std::cerr<<"Start print singleton ... "<<std::endl;
    // print singleton
    for( const auto & pair : reads1 )
        if( common.find(pair.first) == common.end() )
            print_stlfr_reads(pair.second , *o_rs);
    for( const auto & pair : reads2 )
        if( common.find(pair.first) == common.end() )
            print_stlfr_reads(pair.second , *o_rs);
    delete o_rs ;
    std::cerr<<"Done ! "<<std::endl;
    return 0 ;

}
