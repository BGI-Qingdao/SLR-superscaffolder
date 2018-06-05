#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

#include "soap2/fileName.h"
#include "soap2/soap2.h"

#include "stLFR/CBB.h"

#include <vector>
#include <map>

struct AppConfig
{
    struct ConfigBarcodeInfo
    {
        int length;

        std::map<unsigned int , std::vector<unsigned int > > barcodesOnPos;

        std::string format(unsigned int id) const 
        {
            std::ostringstream ost;
            ost<<id<<':'<<length<<'\t';
            for( const auto & i : barcodesOnPos)
            {
                ost<<i.first<<':'<<i.second[0];
                for( int m = 1 ; m < (int)i.second.size(); m++)
                {
                    ost<<'|'<<i.second[m];
                }
                ost<<'\t';
            }
            return ost.str();
        }
    };

    typedef std::map<unsigned int, ConfigBarcodeInfo>  BarcodeOnContig;

    BGIQD::LOG::logger log;

    int bin_size ;

    BarcodeOnContig boc;

    BGIQD::SOAP2::FileNames fName;

    BGIQD::stLFR::BarcodeOnBinArray b2b_array;

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(log,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            auto items = BGIQD::STRING::split(line,"\t");
            boc[std::stoul(items[0])].length = std::stoul(items[1]);
        }
        delete in ;
    }

    void Init(const std::string & prefix,const  int bin)
    {
        bin_size = bin;
        BGIQD::LOG::logfilter::singleton().get("ChopBin",BGIQD::LOG::loglevel::INFO,log);
        fName.Init(prefix);
        b2b_array.Init(1024);
        log<<BGIQD::LOG::lstart()<<"Init finsish ..."<<BGIQD::LOG::lend();
    }

    void LoadBarcodeOnContig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig());
        BGIQD::LOG::timer t(log,"parse all input");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            long readId ;
            unsigned int contigId, pos , barcode;
            char dir;
            std::istringstream ist(line);
            ist>>readId>>contigId>>pos>>dir>>barcode;
            if( boc.find( contigId ) == boc.end() )
            {
                continue;
            }
            boc[contigId].barcodesOnPos[pos].push_back(barcode);
        }
        delete in;
    }

    void PrintBarcodeOnContig()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.BarcodeOnContig());
        if( out == NULL )
            FATAL( "open .barcodeOnContig file to write failed" );

        for( const auto & i : boc )
        {
            (*out)<<i.second.format(i.first)<<std::endl;
        }
        delete out;
    }

    void ChopBin()
    {
        for(const auto & contig :boc )
        {
            unsigned int contigId = contig.first;
            std::map<int,BGIQD::stLFR::BarcodeOnBin> b2b;
            for( const auto & posData : contig.second.barcodesOnPos )
            {
                int pos = posData.first ;
                int binId = pos/bin_size;
                auto & d = b2b[binId];
                d.contigId = contigId;
                d.binId = binId ;

                for( const auto & vv : posData.second)
                {
                    d.collections.IncreaseElement(vv,1);
                }
            }

            for( const auto & i : b2b )
            {
                b2b_array.push_back(i.second);
            }
        }
    }

    void PrintBinInfo()
    {
        BGIQD::stLFR::PrintBarcodeOnBinArray(fName.BarcodeOnBin(),b2b_array);
    }

}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(int , bin_size, "bin size . must be smaller than seed min length");
    DEFINE_ARG_REQUIRED(std::string ,prefix, "prefix . Input xxx.seeds && xxx.read2contig ; Output xxx.barcodeOnBin && xxx.barcodeOnContig");
    DEFINE_ARG_OPTIONAL(bool ,p_b2c , "print barcode on contig ?", "0");
    END_PARSE_ARGS

    config.Init( prefix.to_string() , bin_size.to_int() );

    BGIQD::LOG::timer t(config.log,"ChopBin");

    config.LoadSeeds();

    config.LoadBarcodeOnContig();

    if( p_b2c.to_bool() )
        config.PrintBarcodeOnContig();

    config.ChopBin();

    config.PrintBinInfo();

    return 0;
}
