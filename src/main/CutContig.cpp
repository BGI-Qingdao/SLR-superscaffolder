#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"
#include "common/freq/freq.h"

#include "biocommon/fasta/fasta.h"

struct ContigCutInfo
{
    unsigned int contig ;
    int start ;
    int end ;
};

std::map<unsigned int , ContigCutInfo> data_cache ;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, contig_cut,"the contig_cut file ");
    END_PARSE_ARGS;

    std::string line ;

    auto in = BGIQD::FILES::FileReaderFactory::
        GenerateReaderFromFileName(contig_cut.to_string());
    if( in == NULL )
        FATAL( " failed to read contig_cut file !!!" );
    while( ! std::getline( *in , line).eof() )
    {
        ContigCutInfo info ;
        std::istringstream ist(line) ;
        ist>>info.contig>>info.start>>info.end ;
        data_cache[info.contig] = info ;
    }
    delete in ;

    typedef BGIQD::FASTA::SOAP2ContigHead Header ;
    typedef BGIQD::FASTA::Fasta<Header> Contig ;
    typedef BGIQD::FASTA::FastaReader<Contig> ContigReader ;

    Contig tmp ;
    ContigReader reader ;

    while( reader.LoadNextFasta(std::cin , tmp ) )
    {
        unsigned int cId =  tmp.head.contigId ;
        if( data_cache.find( cId ) == data_cache.end () )
            continue ;
        auto item = data_cache.find( cId )->second;
        if( item.end - item.start + 1  != tmp.head.len )
        {
            tmp.head.len = item.end - item.start + 1  ;
            tmp.seq.atcgs = tmp.seq.atcgs.substr( item.start -1 , tmp.head.len );
        }
        std::cout<<tmp.head.Head()<<'\n';
        std::cout<<tmp.seq.Seq(100);
    }
    return 0;
}
