#include "algorithm/interval/Interval.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

#include "soap2/contigIndex.h"
#include "soap2/fileName.h"
#include "soap2/soap2.h"

#include "stLFR/CBB.h"

#include <vector>
#include <map>
#include <set>
#include <tuple>

struct AppConfig
{
    typedef std::map<unsigned int, BGIQD::stLFR::ContigBarcodeInfo>  BarcodeOnContig;
    typedef BGIQD::INTERVAL::Interval<int,
                BGIQD::INTERVAL::Left_Close_Right_Close > BinInterval;

    BGIQD::LOG::logger log;

    int bin_size ;

    int del_at_tail;

    float bin_factor ;

    bool flatten ;

    enum WorkingMode
    {
        Unknow = 0 ,
        EqualBin = 1 ,
        Head_Tail = 2 ,
        OneBin = 3 ,
        OverlapBin = 4
    };

    WorkingMode work_mode ;

    bool head_tail_only;

    BarcodeOnContig boc;

    std::map<unsigned int ,BGIQD::SOAP2::ContigIndex> seeds;

    BGIQD::SOAP2::FileNames fName;

    BGIQD::stLFR::BarcodeOnBinArray b2b_array;

    std::tuple<bool , int > IsBarcodeInBin( const BinInterval & bin , int barcode_start_pos )
    {
        if( ! flatten )
        {
            bool in = bin.IsContain(barcode_start_pos);
            int len = in ? 1 : 0 ;
            return std::make_tuple( in , len );
        }

        if( barcode_start_pos >= bin.min && barcode_start_pos <= bin.max )
        {
            int total_in_bin = ( bin.max - barcode_start_pos + 1 ) > 100 ?
                100 : ( bin.max - barcode_start_pos + 1 ) ;
            return std::make_tuple(  true , total_in_bin );
        }
        else if ( barcode_start_pos <  bin.min && barcode_start_pos + 99 >= bin.min )
        {
            int right = barcode_start_pos + 99 ;
            if (right > bin.max ) 
                right = bin.max ;
            return std::make_tuple( true , right - bin.min +1 );
        }
        else
        {
            return std::make_tuple( false , 0 );
        }
    }

    int overlap_size  ;
    std::map<int ,BinInterval> MakeBin(int contig_len)
    {
        std::map<int , BinInterval > ret ;
        int contig_used_len = contig_len - del_at_tail;

        if( work_mode == WorkingMode::OneBin)
        {
            ret[0] = BinInterval(-100,contig_used_len+100);
            return  ret ;
        }
        else if ( work_mode == WorkingMode::OverlapBin ) {
            for( int i = 0 ; ( overlap_size * i + bin_size ) <= contig_len ; i++ ) {
                ret[i] = BinInterval(overlap_size*i , overlap_size*i+bin_size-1);
            }
            return  ret ;
        }
        else if ( work_mode == WorkingMode::Head_Tail)
        {
            int half = contig_used_len / 2 ;
            int end1 = half -1 ;
            int start2 = half ;
            if ( end1 > max_bin_size )
            {
                end1 = max_bin_size ;
                start2 = contig_used_len - max_bin_size +1 ;
            }
            if( end1 > bin_size )
                end1 = bin_size ;
            if( start2 < ( contig_used_len - bin_size +1) )
                start2 = contig_used_len - bin_size +1 ;

            ret[0] = BinInterval( 1 , end1 );
            ret[1] = BinInterval( start2 , contig_used_len);
            return ret ;
        }
        else if ( work_mode != WorkingMode::EqualBin )
        {
            assert(0);
        }
        else
        {
            ;
        }

        int start = 1 ;
        int end = contig_used_len ;
        int bin_num = contig_used_len / bin_size + ( (float ( contig_used_len % bin_size ) /(float) bin_size) >= bin_factor? 1 : 0 );
        if ( bin_num == 1 )
        {
            start += (end -start- bin_size) / 2;
            ret[0] = BinInterval(start , start+bin_size -1 );
        }
        else 
        {
            assert( bin_num >= 2 );
            for( int i = 0  ; i< bin_num ; i++ )
            {
                assert( start <= end );
                if( i % 2 == 0 )
                {
                    // Make sure we do not chop bin from very middle area .
                    if ( start > max_bin_size )
                    {
                        break ;
                    }
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
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fName.seeds(middle_name)) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            BGIQD::SOAP2::ContigIndex tmp;
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

        assert( int(work_mode) > 0 );
        assert( int(work_mode) < 4 );
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
                for( const auto & vv : posData.second )
                {
                    for( const auto & pair : bin_intervals )
                    {
                        bool in ; int in_len ;
                        std::tie( in , in_len ) = IsBarcodeInBin( pair.second , pos ) ;
                        //if(pair.second.IsContain(pos))
                        if( in && in_len > 0 )
                        {
                            int binId = pair.first ;
                            auto & d = b2b[binId];
                            d.start = pair.second.min ;
                            d.end = pair.second.max ;
                            d.contigId = contigId;
                            d.binId = binId ;
                            d.collections.IncreaseElement(vv,in_len);
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
        BGIQD::stLFR::PrintBarcodeOnBinArray(fName.BarcodeOnBin(middle_name),b2b_array);
    }

    int max_bin_size ;
    std::string middle_name ;
}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(int , bin_size, "bin size . must be smaller than seed min length");
    DEFINE_ARG_REQUIRED(std::string ,prefix, "prefix . Input xxx.seeds && xxx.barcodeOnContig ; Output xxx.barcodeOnBin");
    //DEFINE_ARG_OPTIONAL(bool ,p_b2c , "print barcode on contig", "0");
    DEFINE_ARG_OPTIONAL(int , work_mode, " the work_mode for chopbin : \n\
                            1 for chop bin with equal bin size \n\
                            2 for chop bin only at contig head and tail \n\
                            3 for chop 1 bin for a contig \n\
                            4 chop bin with certain overlaps "
                            , "1");

    DEFINE_ARG_OPTIONAL(float ,bin_factor , "factor of smallest bin in the middle", "0.5");
    DEFINE_ARG_OPTIONAL(int,  max_bin_size , "max bin area from head & tail " ,"15000");
    DEFINE_ARG_OPTIONAL(bool,  flatten, "flatten mode " ,"false");
    DEFINE_ARG_OPTIONAL(std::string,  middle_name, "the middle name of output suffix " ,"");
    DEFINE_ARG_OPTIONAL(int ,overlap_len, "overlap_len for mode 4", "100");
    END_PARSE_ARGS
    config.overlap_size = overlap_len.to_int();
    config.work_mode = static_cast<AppConfig::WorkingMode>(work_mode.to_int());
    config.max_bin_size = max_bin_size.to_int();
    config.flatten = flatten.to_bool() ;
    config.middle_name = middle_name.to_string() ;
    config.Init( prefix.to_string() , bin_size.to_int() , bin_factor.to_float());

    BGIQD::LOG::timer t(config.log,"ChopBin");

    config.LoadSeeds();

    config.LoadBarcodeOnContig();

    //if( p_b2c.to_bool() )
    //    config.PrintBarcodeOnContig();

    config.ChopBin();

    config.PrintBinInfo();

    return 0;
}
