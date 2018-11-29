#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "biocommon/fasta/fasta.h"
#include "biocommon/fastq/fastq.h"

#include "design_pattern/Publish_Subscribe.h"
#include "design_pattern/Filter.h"

#include "soap2/contigIndex.h"
#include "soap2/fileName.h"

#include <fstream>
#include <map>
#include <vector>
#include <set>

struct AppConfig
{

    typedef BGIQD::FASTA::ScaffSplitGapHead GapHead;
    typedef BGIQD::FASTA::Fasta<GapHead> GapItem;
    typedef BGIQD::FASTA::FastaReader<GapItem> GapReader;

    typedef BGIQD::FASTQ::Fastq<BGIQD::FASTQ::stLFRHeader>  stLFRRead ;
    typedef BGIQD::FASTQ::FastqReader<stLFRRead> stLFRReader ;

    typedef BGIQD::DesignPattern::ISubscriber<stLFRRead> ISub;
    typedef BGIQD::DesignPattern::IPublisher<stLFRRead> IPub;

    typedef BGIQD::DesignPattern::Filter<std::string, IPub> Filter;

    BGIQD::SOAP2::ContigIndexMap contigIndexs;

    Filter read_name_filter;

    Filter barcode_filter ;

    struct ScaffPackage
    {
        int id ;
        std::ofstream * scaff ; 
        std::ofstream * r1; 
        std::ofstream * r2;
        const static int reads_buffer_size = 1024;

        void Init(int id)
        {
            assert(id>0);
            std::string prefix = "tmp_gap_filler_" + std::to_string(id);
            // open 3 files
            scaff = new std::ofstream(prefix+".scaff.fa");
            r1 = new std::ofstream(prefix+".r1.fq");
            r2 = new std::ofstream(prefix+".r2.fq");
            index = 0 ;
            r1_flag = true ;
        }

        void AddScaff(const GapItem & item )
        {
            // save scaffold
            (*scaff)<<item.head.Head()<<'\n';
            (*scaff)<<item.seq.Seq(100);
            delete scaff ;
            scaff = NULL ;
        }
        stLFRRead buffer[reads_buffer_size];
        int index ;

        void CleanBufferR1()
        {
            for( int i = 0 ; i < index ; i ++ )
            {
                (*r1)<<buffer[i].head.Head()<<'\n';
                (*r1)<<buffer[i].seq.Seq(10000000);
            }
            index = 0 ;
        }

        void AddR1( const stLFRRead & read )
        {
            buffer[index] = read ;
            index ++ ;
            if ( index >= reads_buffer_size )
                CleanBufferR1() ;
        }

        void CleanBufferR2()
        {
            for( int i = 0 ; i < index ; i ++ )
            {
                (*r2)<<buffer[i].head.Head()<<'\n';
                (*r2)<<buffer[i].seq.Seq(10000000);
            }
            index = 0 ;
        }

        void AddR2( const stLFRRead & read )
        {
            buffer[index] = read ;
            index ++ ;
            if ( index >= 10000 )
                CleanBufferR2() ;
        }

        void AddRead(  const stLFRRead & read )
        {
            if( r1_flag )
                AddR1( read );
            else
                AddR2( read );
        }

        ~ScaffPackage()
        {
            if ( scaff ) delete  scaff ;
            if ( r1 ) delete r1;
            if ( r2 ) delete r2;
        }

        void EndR1()
        {
            CleanBufferR1();
            r1_flag = false ;
        }

        void EndR2()
        {
            CleanBufferR2();
        }
        bool r1_flag ;
    };

    struct PEGap : public ScaffPackage , public ISub
    {
        BGIQD::SOAP2::ContigIndex c1 ;

        BGIQD::SOAP2::ContigIndex c2 ;

        virtual void update_msg( const stLFRRead & msg ) final
        {
            AddRead(msg);
        }
    };

    PEGap *pe_gaps ;

    struct TrunkGap : public ScaffPackage , public ISub
    {
        BGIQD::SOAP2::ContigIndex c1 ;
        BGIQD::SOAP2::ContigIndex c2 ;

        std::set<std::string> barcodes ;

        std::set<std::string> maped_reads_str;

        std::vector<stLFRRead> maped_reads ;
        std::vector<stLFRRead> unmaped_reads ;
        virtual void update_msg( const stLFRRead & msg ) final
        {
            // PE not barcode keep normal ;
            if( barcodes.find(msg.head.barcode_str) ==  barcodes.end() )
                AddRead(msg);
            // barcode reads save first.
            else
            {
                if( maped_reads_str.find( msg.head.readName ) == maped_reads_str.end())
                    maped_reads.push_back(msg);
                else
                    unmaped_reads.push_back(msg);
            }
        }

        // random fake maped reads & unmaped reads into PE reads
        void EndAllReads()
        {
            for(size_t i = 0 ; i < unmaped_reads.size() ; i++ )
            {
                AddR1(maped_reads[i%maped_reads.size()] );
            }
            CleanBufferR1();
            for(size_t i = 0 ; i < unmaped_reads.size() ; i++ )
            {
                AddR2(unmaped_reads[i] );
            }
            CleanBufferR2();
        }

    };

    TrunkGap * trunk_gaps;

    BGIQD::SOAP2::FileNames fNames;
    BGIQD::LOG::logger loger;

    void Init(const std::string & prefix)
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("Reads2Gap",BGIQD::LOG::loglevel::INFO, loger);
    }

    std::vector<GapItem> gap_buffer ;

    void LoadGaps()
    {

        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.scaff_gap2filler_seqs());
        if( in == NULL )
            FATAL(" failed to open xxx.scaff_gap2filler_seqs for read!!! ");
        GapReader reader ;
        reader.LoadAllFasta(*in , gap_buffer);
        delete in;
    }

    int pe_gap_num ;

    int trunk_gap_num ;

    std::map<unsigned int , std::vector<int> > pe_c2g_index ;

    std::map<unsigned int , std::vector<int> > trunk_c2g_index ;

    void BuildC2GIndex()
    {
        pe_gap_num = 0 ;

        trunk_gap_num = 0 ;

        for( const auto & i : gap_buffer )
        {
            if( i.head.gap_type != GapHead::GapType::TRUNK )
                pe_gap_num ++ ;
            else
                trunk_gap_num ++ ;
        }

        pe_gaps = new PEGap[pe_gap_num] ;
        trunk_gaps = new TrunkGap[trunk_gap_num];

    }

    void LoadContigIndex()
    {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fNames.ContigIndex());
        if(in == NULL)
            FATAL(" failed to open xxx.contigIndex for read!!! ");
        contigIndexs.LoadContigIndexs(*in);
        delete in ;
        contigIndexs.BuildReverseCompleteContigs();
    }

    void BuildGaps()
    {

    }

}config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.minTree ");
    END_PARSE_ARGS

    config.LoadContigIndex();
    config.LoadGaps();
    config.BuildGaps();
    config.LoadPESingles();
    config.LoadBarcodeOnGaps();
    config.ParseReads();
}
