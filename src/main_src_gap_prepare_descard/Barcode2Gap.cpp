#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "algorithm/collection/collection.h"

#include "biocommon/fasta/fasta.h"

#include "soap2/fileName.h"
#include "soap2/contigIndex.h"

#include "stLFR/TrunkGap.h"
#include "stLFR/CBB.h"

#include <set>

struct AppConfig
{
    typedef BGIQD::FASTA::ScaffSplitGapHead Header;

    typedef BGIQD::FASTA::Fasta<Header>  GapFasta;

    typedef BGIQD::Collection::Collection<int> BarcodeCollection;

    BGIQD::SOAP2::FileNames fNames;

    BGIQD::LOG::logger loger;

    std::map<unsigned int , BarcodeCollection > contig_barcode_info;

    std::vector<Header> gap_infos;

    std::set<unsigned int> used_contigs ;

    std::map<int,std::string> barcodes;

    BGIQD::SOAP2::ContigIndexMap contigIndexs ;

    void Init(const std::string & prefix)
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("Barcode2Gap",BGIQD::LOG::loglevel::INFO, loger);
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

    void  LoadScaffoldsGaps()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.scaff_gap2filler_seqs());
        if( in == NULL )
            FATAL(" failed to open xxx.scaff_gap2filler_seqs for read!!! ");
        BGIQD::FASTA::FastaReader<GapFasta> Reader ;

        GapFasta tmp ;

        while(Reader.LoadNextFasta(*in,tmp))
        {
            gap_infos.push_back(tmp.head);
            used_contigs.insert(tmp.head.next_base_contig);
            used_contigs.insert(tmp.head.prev_base_contig);
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

    void PrintBarcode_Intersection_OnGaps()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.barcodeOnGaps_intersection());
        if( out == NULL )
            FATAL(" failed to open xxx.barcodeOnGaps_intersection for write!!! ");

        for( const auto & gap: gap_infos)
        {
            (*out)<<gap.Head()<<'\n';
            auto barcode_both = 
                BarcodeCollection::Intersection(
                        contig_barcode_info[ gap.prev_base_contig]
                    ,   contig_barcode_info[ gap.next_base_contig] );

            for(const auto & pair : barcode_both.elements)
            {
                (*out)<<barcodes.at(pair.first)<<'\n';
            }
        }
        delete out;
    }

    void PrintBarcode_Union_OnGaps()
    {
        auto out = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.barcodeOnGaps_union());
        if( out == NULL )
            FATAL(" failed to open xxx.barcodeOnGaps_union for write!!! ");

        for( const auto & gap: gap_infos)
        {
            (*out)<<gap.Head()<<'\n';
            auto barcode_both = 
                BarcodeCollection::Union(
                        contig_barcode_info[ gap.prev_base_contig]
                    ,   contig_barcode_info[ gap.next_base_contig] );

            for(const auto & pair : barcode_both.elements)
            {
                (*out)<<barcodes.at(pair.first)<<'\n';
            }
        }
        delete out;
    }

    void LoadBarcodeOnContig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fNames.BarcodeOnContig());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnContig to read !");
        BGIQD::LOG::timer t(loger,"LoadBarcodeOnContig");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            if( used_contigs.find( tmp.contig_id) == used_contigs.end() )
            {
                continue;
            }
            BarcodeCollection tdata;
            for(const auto & pair : tmp.barcodesOnPos)
            {
                const auto & barcodes = pair.second ;
                for(auto b : barcodes)
                {
                    tdata.IncreaseElement(b);
                }
            }
            contig_barcode_info[tmp.contig_id] = tdata;
        }
        delete in;
    }

}config;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix  ");
    END_PARSE_ARGS

    config.Init(prefix.to_string());
    config.LoadContigIndex();
    config.LoadBarcodeList();
    config.LoadScaffoldsGaps();
    config.LoadBarcodeOnContig();
    config.PrintBarcode_Intersection_OnGaps();
    config.PrintBarcode_Union_OnGaps();
}
