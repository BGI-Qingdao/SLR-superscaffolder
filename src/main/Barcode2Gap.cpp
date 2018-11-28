#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"


#include "biocommon/fasta/fasta.h"

#include "soap2/fileName.h"
#include "stLFR/TrunkGap.h"
#include "stLFR/CBB.h"

#include <set>

struct AppConfig
{
    BGIQD::SOAP2::FileNames fNames;
    BGIQD::LOG::logger loger;

    typedef BGIQD::stLFR::TrunkGap<int> GapInfo;

    typedef BGIQD::FASTA::ScaffSplitGapHead Header;

    typedef BGIQD::FASTA::Fasta<Header>  GapFasta;

    std::map<int,std::vector<GapInfo>> gaps;
    std::map<int,BGIQD::stLFR::ContigOnBarcode>  c2bs;

    //std::map<int,std::set<int>> scaff2contig;
    std::map<int,std::set<int> > contigOnScaffld;
    typedef BGIQD::FREQ::Freq<std::string> BarcodesInfo;
    std::map<int,BarcodesInfo> barcodeOnscaff;
    std::map<int,std::string> barcodes;

    void Init(const std::string & prefix)
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadContigOnBarcode()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.contigOnBarcode());
        if( in == NULL )
            FATAL(" failed to open xxx.contigOnBarcode for read!!! ");
        auto read_line = [this](const std::string & line)
        {
            BGIQD::stLFR::ContigOnBarcode tmp ;
            tmp.InitFromString(line);
            auto barcode_str = barcodes[tmp.barcode_id];
            for(auto pair : tmp.contig_data)
            {
                int contigId = pair.first;
                if( contigOnScaffld.find( contigId ) == contigOnScaffld.end())
                    continue ;
                for( int scaffoldId : contigOnScaffld[contigId])
                {
                    barcodeOnscaff[scaffoldId].Touch(barcode_str , pair.second );
                }
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,read_line);
        delete in ;
    }

    void  LoadScaffoldsGaps()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.scaff_gap2filler_seqs());
        if( in == NULL )
            FATAL(" failed to open xxx.scaff_gap2filler_seqs for read!!! ");
        BGIQD::FASTA::FastaReader<GapFasta> Reader ;

        GapFasta tmp ;
        int scaffold_id = 1 ;
        while(Reader.LoadNextFasta(*in,tmp))
        {
            if(tmp.head.gap_type != Header::GapType::TRUNK )
                continue ;
            contigOnScaffld[tmp.head.prev_contig].insert(scaffold_id);
            contigOnScaffld[tmp.head.next_contig].insert(scaffold_id);
            scaffold_id ++ ;
        }

        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load gaps done "<<BGIQD::LOG::lend() ;
    }

    void LoadBarcodeList()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.barcodeList());
        if( in == NULL )
            FATAL(" failed to open xxx.barcodeList for read !!! ");
        auto read_line = [this](const std::string & line)
        {
            auto items = BGIQD::STRING::split(line,'\t');
            assert(items.size() == 2 );
            barcodes[std::stoi(items[1])] = items[0];
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,read_line);
        delete in ;
    }

    void PrintBarcodeOnScaffold()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.barcodeOnScaff());
        if( out == NULL )
            FATAL(" failed to open xxx.barcodeOnScaffold for write!!! ");
        for( const auto & pair : barcodeOnscaff )
        {
            (*out)<<">"<<pair.first<<'\n';
            (*out)<<pair.second.ToString();
        }
        delete out;
    }
}config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix , Input xxx.cluster . Output xxx.minTree ");
    END_PARSE_ARGS

    config.Init(prefix.to_string());

    config.LoadBarcodeList();
    config.LoadScaffoldsGaps();
    config.LoadContigOnBarcode();
    config.PrintBarcodeOnScaffold();
}
