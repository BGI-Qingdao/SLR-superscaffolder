#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"

#include "biocommon/fasta/fasta.h"
#include "biocommon/seq/tool_func.h"

#include "soap2/fileName.h"

#include "stLFR/ScaffInfo.h"

#include <map>
#include <set>

struct AppConfig
{
    BGIQD::SOAP2::FileNames fNames;

    BGIQD::LOG::logger loger;

    BGIQD::stLFR::ScaffInfoHelper scaff_helper;

    typedef BGIQD::FASTA::SOAP2ContigHead Header;

    typedef BGIQD::FASTA::Fasta<Header> ContigFasta;

    typedef BGIQD::FASTA::FastaReader<ContigFasta> Reader;

    std::map<unsigned int , ContigFasta> contigs;

    std::set<unsigned int> used;

    void Init( const std::string & prefix )
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("ScaffInfo2Seqs",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadScaffInfos()
    {
        auto  in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.scaff_infos());
        if( in == NULL )
            FATAL("failed to open xxx.scaff_infos to read");

        scaff_helper.LoadAllScaff(*in);
        delete in ;
    }

    void LoadAllFasta()
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.contig());
        if( in == NULL )
            FATAL("failed to open xxx.scaff_infos to read");

        Reader reader;
        ContigFasta tmp ;
        while( reader.LoadNextFasta(*in , tmp) )
        {
            contigs[tmp.head.contigId] = tmp ;
        }
        delete in ;
    }

    void PrintScaffGapSeqs()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.scaff_gap2filler_seqs()) ;

        auto get_atcg_withoutn = [&] ( const BGIQD::stLFR::ContigDetail & detail )
        {
            std::string str = contigs.at(detail.contig_id).seq.atcgs ;
            if( ! detail.orientation )
                str = BGIQD::SEQ::seqCompleteReverse(str);
            if( detail.gap_size > 0 )
            { ;}
            else if ( detail.gap_size < 0 )
            {
                assert( (int)str.length() > (int)std::abs(detail.gap_size) );
                str = str.substr(0,str.length() + detail.gap_size);
            }
            else
            { ; }

            return str ;

        };

        auto get_atcg = [&] ( const BGIQD::stLFR::ContigDetail & detail )
        {
            std::string str = contigs.at(detail.contig_id).seq.atcgs ;
            if( ! detail.orientation )
                str = BGIQD::SEQ::seqCompleteReverse(str);
            if( detail.gap_size > 0 )
                str += std::string(detail.gap_size,'N');
            else if ( detail.gap_size < 0 )
            {
                assert( (int)str.length() > (int)std::abs(detail.gap_size) );
                str = str.substr(0,str.length() + detail.gap_size);
                str += "N" ;
            }
            else
            { str += "N" ; }

            return str ;
        };

        for( const auto & pair : scaff_helper.all_scaff)
        {
            const auto & scaff_vec = pair.second.a_scaff ;
            for( int i = 1 ; i < (int)scaff_vec.size() ; i++ )
            {
                const auto & prev = scaff_vec.at(i-1);
                const auto & next = scaff_vec.at(i);
                BGIQD::FASTA::ScaffSplitGapHead header ;
                if( prev.gap_size < 50 )
                {
                    header.gap_type = BGIQD::FASTA::ScaffSplitGapHead::GapType::PE_TRUNK ;
                }
                else
                {
                    header.gap_type = BGIQD::FASTA::ScaffSplitGapHead::GapType::TRUNK ;
                }
                header.gap_index = i ;
                header.next_base_contig = next.contig_id ;
                header.prev_base_contig = prev.contig_id ;
                header.next_contig = next.orientation ? next.contig_id : next.contig_id +1 ;
                header.prev_contig = prev.orientation ? prev.contig_id : prev.contig_id +1 ;
                header.scaff_id = pair.first ;
                (*out)<<header.Head()<<'\n';
                std::string str ;
                str +=get_atcg(prev);
                str += get_atcg_withoutn(next);
                (*out)<<BGIQD::SEQ::blockSeq(str,100);
            }
        }
        delete out ;
    }

} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.scaff_infos; Output xxx.scaff_seqs");
    END_PARSE_ARGS;

    config.Init( prefix.to_string() );
    BGIQD::LOG::timer t(config.loger,"ScaffInfo2Seqs");

    config.LoadScaffInfos();
    config.LoadAllFasta();
    config.PrintScaffGapSeqs();
    return 0;
}
