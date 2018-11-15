#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "biocommon/fasta/fasta.h"
#include "biocommon/fastq/fastq.h"

#include "soap2/fileName.h"
#include "stLFR/TrunkGap.h"
#include "stLFR/CBB.h"

#include <fstream>

struct AppConfig
{

    BGIQD::SOAP2::FileNames fNames;
    BGIQD::LOG::logger loger;

    std::string r1;
    std::string r2;

    typedef BGIQD::FASTA::Fasta<BGIQD::FASTA::NormalHead> ScaffItem;
    typedef BGIQD::FASTA::FastaReader<ScaffItem> Reader;
    std::vector<ScaffItem> scaff_buffer;

    typedef BGIQD::FREQ::Freq<std::string> BarcodesInfo;
    std::map<int,BarcodesInfo> barcodeOnscaff ;

    typedef BGIQD::FASTQ::Fastq<BGIQD::FASTQ::stLFRHeader>  stLFRRead ;
    typedef BGIQD::FASTQ::FastqReader<stLFRRead> stLFRReader ;

    bool IsBarcodeInScaff( const std::pair<int,BarcodesInfo> & info , const std::string & barcode)
    {
        return info.second.data.find(barcode) != info.second.data.end() ;
    }

    std::set<int> FilterScaffolds( const std::string & barcode )
    {
        std::set<int> ret ;
        for( const auto & pair : barcodeOnscaff )
        {
            if( IsBarcodeInScaff( pair,barcode ) )
                ret.insert( pair.first );
        }
        return ret ;
    }

    void Init(const std::string & prefix
            , const std::string & read1
            , const std::string & read2 
            )
    {
        r1 = read1 ;
        r2 = read2 ;
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadBarcodeConScaffold()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.barcodeOnScaff());
        if( in == NULL )
            FATAL(" failed to open xxx.barcodeOnscaf for read!!! ");

        int curr_scaff = -1 ;
        auto read_line = [this, &curr_scaff](const std::string & line)
        {
            if( line[0] == '>' )
                curr_scaff = std::stoi(line);
            else
            {
                auto items = BGIQD::STRING::split(line,'\t');
                assert(items.size() == 2 );
                assert(curr_scaff != -1 );
                barcodeOnscaff[curr_scaff].Touch(items[0],std::stoi(items[1]));
            }
        };

        BGIQD::FILES::FileReaderFactory::EachLine(*in,read_line);
        delete in ;
    }

    void LoadScaff_Seqs()
    {
        auto in  = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.scaff_seqs());
        if( in == NULL )
            FATAL(" failed to open xxx.scaff_seqs for read!!! ");
        Reader::LoadAllFasta(*in , scaff_buffer);
        delete in ;
    }


    struct ScaffPackage
    {
        int id ;
        std::ofstream * scaff ; 
        std::ofstream * r1; 
        std::ofstream * r2;

        void Init(int id)
        {
            assert(id>0);
            std::string prefix = "tmp_gap_filler_" + std::to_string(id);
            // open 3 files
            scaff = new std::ofstream(prefix+".scaff.fa");
            r1 = new std::ofstream(prefix+".r1.fq");
            r2 = new std::ofstream(prefix+".r2.fq");
        }

        void AddScaff(const ScaffItem & item )
        {
            // save scaffold
            (*scaff)<<item.head.Head()<<'\n';
            (*scaff)<<item.seq.Seq(100);
            delete scaff ;
            scaff = NULL ;
        }

        void AddR1( const stLFRRead & read )
        {
            (*r1)<<read.head.Head()<<'\n';
            (*r1)<<read.seq.Seq(10000000);
        }

        void AddR2( const stLFRRead & read )
        {
            (*r2)<<read.head.Head()<<'\n';
            (*r2)<<read.seq.Seq(10000000);
        }

        ~ScaffPackage()
        {
            if ( scaff ) delete  scaff ;
            if ( r1 ) delete r1;
            if ( r2 ) delete r2;
        }

    };
    std::map<int , ScaffPackage> scaffPrinters;

    void PrepareFiles()
    {
        for( const auto & pair : barcodeOnscaff )
        {
            scaffPrinters[pair.first].Init(pair.first);
            scaffPrinters[pair.first].AddScaff(scaff_buffer[pair.first-1]);
        }
    }

    void ParseRead1()
    {
        auto in  = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(r1);
        if( in == NULL )
            FATAL(" failed to open r1 for read!!! ");

        stLFRReader reader;
        stLFRRead tmp ;
        while ( reader.LoadNextFastq(*in , tmp) )
        {
            auto rets = FilterScaffolds(tmp.head.barcode_str);

            for( const auto & item : rets )
            {
                scaffPrinters[item].AddR1(tmp);
            }
        }
        delete in ;
    }

    void ParseRead2()
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(r2);
        if( in == NULL )
            FATAL(" failed to open r1 for read!!! ");

        stLFRReader reader ;
        stLFRRead tmp ;
        while ( reader.LoadNextFastq(*in , tmp) )
        {
            auto rets = FilterScaffolds(tmp.head.barcode_str);

            for( const auto & item : rets )
            {
                scaffPrinters[item].AddR2(tmp);
            }
        }
        delete in ;
    }

    void ParseAllReads()
    {
        ParseRead1();
        ParseRead2();
    }

} config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.minTree ");
        DEFINE_ARG_REQUIRED(std::string , r1, "read1");
        DEFINE_ARG_REQUIRED(std::string , r2, "read2");
    END_PARSE_ARGS

    config.Init(prefix.to_string(), r1.to_string() , r2.to_string());

    config.LoadBarcodeConScaffold();

    config.LoadScaff_Seqs();

    config.PrepareFiles();

    config.ParseAllReads();
}
