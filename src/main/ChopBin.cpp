#include "algorithm/interval/Interval.h"

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
#include <set>

struct AppConfig
{
    typedef std::map<unsigned int, BGIQD::stLFR::ContigBarcodeInfo>  BarcodeOnContig;
    typedef BGIQD::INTERVAL::Interval<int,
                BGIQD::INTERVAL::Left_Close_Right_Close > BinInterval;

    BGIQD::LOG::logger log;

    int bin_size ;

    int del_at_tail;

    float bin_factor ;

    bool head_tail_only;

    BarcodeOnContig boc;

    std::map<unsigned int ,BGIQD::stLFR::ContigIndex> seeds;

    BGIQD::SOAP2::FileNames fName;

    BGIQD::stLFR::BarcodeOnBinArray b2b_array;

    std::map<int ,BinInterval> MakeBin(int contig_len)
    {
        std::map<int , BinInterval > ret ;
        int contig_used_len = contig_len - del_at_tail;
        int bin_num = contig_used_len / bin_size + ( (float ( contig_used_len % bin_size ) /(float) bin_size) >= bin_factor? 1 : 0 );
        int start = 1 ;
        int end = contig_used_len ;

        if ( bin_num == 1 )
        {
            start += (end -start- bin_size) / 2;
            ret[0] = BinInterval(start , start+bin_size -1 );
        }
        else 
        {
            assert( bin_num >= 2 );
            if( head_tail_only )
                bin_num = 2;
            for( int i = 0  ; i< bin_num ; i++ )
            {
                assert( start <= end );
                if( i % 2 == 0 )
                {
                    ret[i/2] = 
                        BinInterval(start , start+bin_size -1 );
                    start += bin_size ;
                }
                else
                {
                    ret[bin_num-(i/2)-1] = 
                        BinInterval(end -bin_size +1 , end );
                    end -= bin_size ;
                }
            }
        }
        assert( bin_num > 0 );
        return ret ;
    };


    void LoadSeeds()
    {
        BGIQD::LOG::timer t(log,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            BGIQD::stLFR::ContigIndex tmp;
            tmp.InitFromString(line);
            seeds[tmp.contig] = tmp;
            boc[tmp.contig].contig_id= tmp.contig;
        }
        delete in ;
    }

    void Init(const std::string & prefix
            ,const  int bin 
            , float bin_f)
    {
        bin_size = bin;
        del_at_tail = 0;
        bin_factor = bin_f;
        BGIQD::LOG::logfilter::singleton().get("ChopBin",BGIQD::LOG::loglevel::INFO,log);
        fName.Init(prefix);
        b2b_array.Init(1024);
        log<<BGIQD::LOG::lstart()<<"Init finsish ..."<<BGIQD::LOG::lend();
    }

    void LoadBarcodeOnContig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.BarcodeOnContig());
        if(! in )
            FATAL( "failed to open xxx.barcodeOnContig to read !");
        BGIQD::LOG::timer t(log,"parse all input");
        std::string line;
        while(!std::getline(*in,line).eof())
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            if( boc.find( tmp.contig_id) == boc.end() )
            {
                continue;
            }
            boc[tmp.contig_id] = tmp ;
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
            (*out)<<i.second.ToString()<<std::endl;
        }
        delete out;
    }


    void ChopBin()
    {
        for(const auto & contig :boc )
        {
            unsigned int contigId = contig.first;
            int len = seeds[contigId].length ;
            auto bin_intervals = MakeBin(len);
            std::map<int,BGIQD::stLFR::BarcodeOnBin> b2b;
            for( const auto & posData : contig.second.barcodesOnPos )
            {
                int pos = posData.first ;
                for( const auto & vv : posData.second)
                {
                    for( const auto & pair : bin_intervals )
                    {
                        if(pair.second.IsContain(pos))
                        {
                            int binId = pair.first ;
                            auto & d = b2b[binId];
                            d.start = pair.second.min ;
                            d.end = pair.second.max ;
                            d.contigId = contigId;
                            d.binId = binId ;
                            d.collections.IncreaseElement(vv,1);
                        }
                    }
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
    //DEFINE_ARG_REQUIRED(int , delete_tail, "delete size at contig tail . depends on you alignment tools and args");
    DEFINE_ARG_OPTIONAL(float ,bin_factor , "factor of smallest bin in the middle", "0.5");
    DEFINE_ARG_REQUIRED(std::string ,prefix, "prefix . Input xxx.seeds && xxx.barcodeOnContig ; Output xxx.barcodeOnBin");
    //DEFINE_ARG_OPTIONAL(bool ,p_b2c , "print barcode on contig", "0");
    DEFINE_ARG_OPTIONAL(bool ,ht_only, "only chop bin at head and tail", "0");
    END_PARSE_ARGS

    config.Init( prefix.to_string() , bin_size.to_int() , bin_factor.to_float());
    config.head_tail_only= ht_only.to_bool() ;

    BGIQD::LOG::timer t(config.log,"ChopBin");

    config.LoadSeeds();

    config.LoadBarcodeOnContig();

    //if( p_b2c.to_bool() )
    //    config.PrintBarcodeOnContig();

    config.ChopBin();

    config.PrintBinInfo();

    return 0;
}
