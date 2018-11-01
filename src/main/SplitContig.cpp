#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"

#include "soap2/fileName.h"
#include "biocommon/fasta/fasta.h"

struct AppConfig {

    BGIQD::SOAP2::FileNames fNames;

    typedef BGIQD::FASTA::Fasta<BGIQD::FASTA::NormalHead> Fasta;

    typedef BGIQD::FASTA::FastaReader<Fasta> Reader ;

    std::vector<Fasta> contig_buffer;

    int threshold ;

    void Init(const std::string & prefix , int t)
    {
        fNames.Init(prefix);
        threshold = t ;
    }

    void LoadConfig()
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName( fNames.contig() );
        if( in == NULL )
            FATAL( "failed to open xxx.config for read !!! ");
        Reader::LoadAllFasta(*in, contig_buffer );
        delete in ;
    }

    void PrintShort() 
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.contig_short());
        if( out == NULL )
            FATAL( "failed to open xxx.config_short for write!!! ");
        for( const auto & x : contig_buffer)
        {
            if( x.seq.Len() < threshold )
            {
                (*out)<<x.head.Head()<<'\n';
                (*out)<<x.seq.Seq(100);
            }
        }
        delete out ;
    }

    void PrintLong() 
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.contig_long());
        if( out == NULL )
            FATAL( "failed to open xxx.config_long for write !!! ");
        for( const auto & x : contig_buffer)
        {
            if( x.seq.Len() >= threshold )
            {
                (*out)<<x.head.Head()<<'\n';
                (*out)<<x.seq.Seq(100);
            }
        }
        delete out ;
    }
} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.contig ; Output xxx.contig_long && xxx.contig_short");
        DEFINE_ARG_REQUIRED(int , threshold,"smaller than threshold goes to short , otherwise to long");
    END_PARSE_ARGS;

    config.Init( prefix.to_string() , threshold.to_int() );
    config.LoadConfig();
    config.PrintShort();
    config.PrintLong();

    return 0;
}
