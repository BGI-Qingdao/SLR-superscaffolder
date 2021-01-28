#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/misc/Error.h"

#include "utils/seq/fasta.h"
#include "utils/seq/tool_func.h"
#include "utils/agp/agp.h"

#include "utils/misc/fileName.h"

#include "stLFR/ScaffInfo.h"
#include "stLFR/StringIdCache.h"

#include <map>
#include <set>
#include <sstream>

struct AppConfig
{
    BGIQD::MISC::FileNames fNames;

    BGIQD::LOG::logger loger;

    BGIQD::stLFR::ScaffInfoHelper scaff_helper;

    typedef BGIQD::SEQ::SOAP2ContigHead Header;
    typedef BGIQD::SEQ::Fasta<Header> ContigFasta;
    typedef BGIQD::SEQ::FastaReader<ContigFasta> Reader;

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
            FATAL("failed to open xxx.contig to read");

        Reader reader;
        ContigFasta tmp ;
        while( reader.LoadNextFasta(*in , tmp) )
        {
            contigs[tmp.head.contigId] = tmp ;
        }
        delete in ;
    }


    void PrintScaffSeqs()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.scaff_seqs()) ;

        auto get_atcg = [&] ( const BGIQD::stLFR::ContigDetail & detail
                , BGIQD::AGP::Scaff2AGPItem & s2a ) -> std::string
        {
            std::string str = contigs.at(detail.contig_id).seq.atcgs ;
            if( ! detail.orientation )
                str = BGIQD::SEQ::seqCompleteReverse(str);
            if( detail.gap_size > 0 )
            {
                s2a.AddSeq(
                        contig_name_cache.Id(detail.contig_id)
                        , 1
                        , detail.contig_len 
                        , (detail.orientation ? '+' : '-' )
                        ) ;

                if( detail.extra.find(BGIQD::stLFR::ContigDetail::ONT_FILL) 
                        != detail.extra.end() )
                {
                    str += detail.extra.at(BGIQD::stLFR::ContigDetail::ONT_FILL);
                    s2a.AddN(detail.extra.at(BGIQD::stLFR::ContigDetail::ONT_FILL).size());
                }
                else
                {
                    int gap_size = detail.gap_size > min_n ?  detail.gap_size : min_n ;
                    str += std::string(gap_size,'N');
                    s2a.AddN(gap_size);
                }
            }
            else if ( detail.gap_size < 0 )
            {
                if(  (int)str.length() +  (int)std::abs(detail.gap_size) > 0 )
                {
                    str = str.substr(0,str.length() + detail.gap_size);
                    s2a.AddSeq(
                            contig_name_cache.Id(detail.contig_id)
                            , 1
                            , detail.contig_len + detail.gap_size
                            , (detail.orientation ? '+' : '-' )
                            ) ;
                }
                str += std::string(min_c,'N');
                s2a.AddN(min_c);
            }
            else
            {
                // gap size = 0
                s2a.AddSeq(
                        contig_name_cache.Id(detail.contig_id)
                        , 1
                        , detail.contig_len 
                        , (detail.orientation ? '+' : '-' )
                        ) ;
                ;
            }
            return str ;
        };

        for( const auto & pair : scaff_helper.all_scaff)
        {
            BGIQD::AGP::Scaff2AGPItem s2a;
            s2a.InitName("scaffold_" + std::to_string(pair.first));
            (*out)<<">scaffold_"<<pair.first<<'\n';
            std::string str ;
            for( const auto & i : pair.second.a_scaff) 
            {
                used.insert(i.contig_id);
                str+= get_atcg(i,s2a);
            }
            agp_cache.data.insert(agp_cache.data.end() ,s2a.Items().begin() ,s2a.Items().end());
            (*out)<<BGIQD::SEQ::blockSeq(str,100);
        }
        for( const auto & pair : contigs )
        {
            if( used.find( pair.first ) == used.end () )
            {
                BGIQD::AGP::Scaff2AGPItem tmp ;
                tmp.InitName(contig_name_cache.Id(pair.second.head.contigId));
                tmp.AddSeq(contig_name_cache.Id(pair.second.head.contigId),
                        1,
                        pair.second.head.len,
                        '+');
                if( contig_name_cache.HasId(pair.second.head.contigId) )
                    (*out)<<'>'<<contig_name_cache.Id(pair.second.head.contigId) <<'\n';
                else
                    (*out)<<pair.second.head.Head()<<'\n';
                (*out)<<pair.second.seq.Seq(100);
                agp_cache.data.push_back(*tmp.Items().begin());
            }
        }
        delete out ;
    }
    int min_n ;
    int min_c ;

    BGIQD::AGP::AGPFile agp_cache;
    void PrintAGP()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.scaff_agp());
        if( out == NULL )
            FATAL( " failed to open xxx.scaff_agp to write !!! ");
        agp_cache.Print(*out);
        delete out ;
    }

    BGIQD::stLFR::IdStringCache contig_name_cache;
    void LoadFakeSOAPNameMap()
    {
        auto fake = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName( "fakesoap.name2index.map.txt" );
        if( fake == NULL )
            INFO_RETURN(" no fakesoap.name2index.map.txt file !!! ");
        delete fake ;
        contig_name_cache.preload = true ;
        contig_name_cache.LoadStringIdCache("fakesoap.name2index.map.txt");
    }
} config;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.scaff_infos; Output xxx.scaff_seqs");
        DEFINE_ARG_OPTIONAL(int, min_n,"min N size for gap in scaffold if not filled ","11");
        DEFINE_ARG_OPTIONAL(int, min_c,"min N size for contig overlap cut ","11");
    END_PARSE_ARGS;
    config.min_c = min_c.to_int() ;
    config.min_n = min_n.to_int();
    config.Init( prefix.to_string() );
    BGIQD::LOG::timer t(config.loger,"ScaffInfo2Seqs");

    config.LoadFakeSOAPNameMap();
    config.LoadScaffInfos();
    config.LoadAllFasta();
    config.PrintScaffSeqs();
    config.PrintAGP();
    return 0;
}
