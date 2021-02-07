/**********************************************************
 *
 * @Brief  :
 *
 *   Parse the read names of read1 and collect informations:
 *      + map string readnames to numbers
 *      + map string barcodes to numbers
 *      + filter low quality barcodes.
 *
 * *******************************************************/
#include "utils/log/log.h"
#include "utils/args/argsparser.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"
#include "utils/misc/Error.h"
#include "utils/string/stringtools.h"
#include "utils/misc/freq.h"
#include "utils/misc/fileName.h"
#include "utils/misc/TagId.h"
#include "utils/misc/flags.h"

#include "stLFR/stLFRRead.h"

#include <set>
#include <string>
#include <cassert>
#include <vector>
#include <iostream>

// define a stLFR read in fastq format
struct Fastq
{
    typedef BGIQD::stLFR::stLFRHeader Header;

    FLAGS_INT ;

    ADD_A_FLAG(1,UnSet);
    ADD_A_FLAG(2,Set_head);
    ADD_A_FLAG(3,Set_seq);
    ADD_A_FLAG(4,Set_3);
    ADD_A_FLAG(5,Set_quality);

    void Reset() 
    {
        head.Reset();
        flags = 0 ;
        Set_UnSet();
    }

    Header head;

    void AddHead(const std::string & line)
    {
        Clean_UnSet();
        Set_Set_head();
        head.Init(line);
    }

    void AddSeq(const std::string & /*line*/ )
    {
        if( Is_UnSet() || ! Is_Set_head() )
        {
            assert(0);
        }
        Set_Set_seq();
    }

    void Add3(const std::string & )
    {
        if( Is_UnSet()
                || ! Is_Set_head() 
                || ! Is_Set_seq() )
        {
            assert(0);
        }
        Set_Set_3();
    }

    void AddQuality( const std::string & /*line*/ )
    {

        if( Is_UnSet() || ! Is_Set_3() )
        {
            assert(0);
        }
        Set_Set_quality();
    }

    bool Is_Setted() const { 
        return  (!Is_UnSet() )
            && Is_Set_head() 
            && Is_Set_seq()
            && Is_Set_3() 
            && Is_Set_quality();
    }
};
// wrap basic functions for loading and parsing fastq
struct FastqReader
{
    static bool IsHead(const std::string & line)
    {
        return ( ! line.empty()) && line[0] == '@' ;
    }
    static bool Is_3(const std::string & line)
    {
        return ( ! line.empty()) && line[0] == '+' ;
    }

    static void LoadAllFastq( 
            std::istream & ist 
            , std::vector<Fastq> & buffer )
    {
        std::string line ;
        Fastq fq;
        fq.Reset();
        long index = 0 ;
        while ( ! std::getline(ist,line).eof() )
        {
            index ++ ;
            if( index % 4 == 1 )
            {
                assert(IsHead(line) && "fastq format invalid!!!");
                if(fq.Is_Setted())
                {
                    buffer.push_back(fq);
                }
                fq.Reset();
                fq.AddHead(line);
            }
            else if ( index % 4 == 2 )
            {
                fq.AddSeq(line);
            }
            else if ( index % 4 == 3 )
            {
                assert( Is_3(line) && "fastq format invalid!!!" );
                fq.Add3(line);
            }
            else 
            {
                fq.AddQuality(line);
            }
        }
        if(fq.Is_Setted())
        {
            buffer.push_back(fq);
        }
        fq.Reset();
    }

    static bool LoadNextFastq(std::istream & ist , Fastq & fq)
    {
        std::string line ;
        int index = 0;
        fq.Reset();
        while ( ! std::getline(ist,line).eof() )
        {
            index ++ ;
            if( index % 4 == 1 )
            {
                assert(IsHead(line) && "fastq format invalid!!!");
                fq.AddHead(line);
            }
            else if ( index % 4 == 2 )
            {
                fq.AddSeq(line);
            }
            else if ( index % 4 == 3 )
            {
                assert( Is_3(line) && "fastq format invalid!!!" );
                fq.Add3(line);
            }
            else 
            {
                fq.AddQuality(line);
            }
            if( index == 4 )
                break ;
        }
        return fq.Is_Setted();
    }
};

/**********************************************************
 *
 *  Define a conf class to wrap all global variables and 
 *  logical funtions.
 *
 * ********************************************************/
struct AppConfig
{
    typedef BGIQD::stLFR::stLFRHeader::ReadType Type ;
    typedef FastqReader Reader;

    BGIQD::LOG::logger loger;

    BGIQD::MISC::StringIdCache barcode_list;

    BGIQD::MISC::StringIdCache unknow_barcode_list;

    BGIQD::MISC::Freq<std::string>  barcode_freq;

    int min_pair_threshold ;

    int max_pair_threshold ;

    std::string read1;

    BGIQD::MISC::FileNames fNames;

    void Init(const std::string & in , const std::string & prefix)
    {
        read1 = in ;
        fNames.Init(prefix);
        loger.Init("ParseReadName");
    }

    // loading read1 line by line and parse the read names.
    void ParseRead1()
    {
        long long index = 1 ;
        auto out = BGIQD::FILES::
            FileWriterFactory::GenerateWriterFromFileName
            (fNames.readNameList());

        BGIQD::LOG::timer timer(loger,"ParseRead1");

        if( out == NULL )
            FATAL( " failed to open xxx.readNameList for write !!! ");

        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(read1);
        if ( in == NULL )
            FATAL( " failed to open read1 for read !!! ");
        Fastq data;
        long next_id =1 ;
        std::set<std::string> barcodes;
        barcodes.insert("0");
        barcodes.insert("0_0");
        barcodes.insert("0_0_0");
        barcodes.insert("0_0_0_0");
        while( Reader::LoadNextFastq(*in , data))
        {
            (*out)<<data.head.readName<<'\t'<<index<<'\n';
            index ++ ;
            if( index >= 1000000 && index % 1000000 == 1 )
            {
                loger<<BGIQD::LOG::lstart()<<"process "
                    <<index<<" read... "<<BGIQD::LOG::lend();
            }

            if( data.head.type ==
                    Type::readName_barcodeStr_index_barcodeNum
            ||  data.head.type ==  Type::readName_barcodeStr
            ||  data.head.type == Type::readName_barcodeStr_index
                    )

            {
                barcode_list.preload = true ;
                if( barcodes.find(data.head.barcode_str)
                 == barcodes.end() )
                {
                    barcode_list.data.AssignTag(
                            data.head.barcode_str,
                            next_id);
                    next_id ++ ;
                    barcodes.insert(data.head.barcode_str);
                }
            }
            else
            {
                unknow_barcode_list.preload = true ;
                unknow_barcode_list.Id(data.head.barcode_str);
            }
            barcode_freq.Touch(data.head.barcode_str);
        }
        delete in ;
        delete out ;
    }

    // print barcode informations.
    void PrintBarcodeList()
    {
        if( barcode_list.preload == unknow_barcode_list.preload )
        {
            loger<<BGIQD::LOG::lstart()
                <<"ERROR :  some has barcode "<<barcode_list.preload
                <<" some has no barcode "<<unknow_barcode_list.preload
                <<BGIQD::LOG::lend();
        }
        BGIQD::LOG::timer timer(loger,"PrintBarcodeList");
        if( barcode_list.preload )
        {
            loger<<BGIQD::LOG::lstart()<<"has barcode num in read1"
                <<BGIQD::LOG::lend();
            barcode_list.Print(fNames.barcodeList());
        }
        else
        {
            loger<<BGIQD::LOG::lstart()<<"no barcode num in read1"
                <<BGIQD::LOG::lend();
            unknow_barcode_list.Print(fNames.barcodeList());
        }
    }

    // mask low quality barcode as zero-barcode.
    void MaskTooLowBarcode()
    {
        for( const auto & pair : barcode_freq.data )
        {
            long num = pair.second ;
            const std::string & barcode_str = pair.first ;
            if( num < min_pair_threshold )
                barcode_list.data.AssignTag(barcode_str,0);
        }
    }

    // mask too density barcode as zero-barcode because they may caused by barcode collision.
    void MaskTooHighBarcode()
    {
        for( const auto & pair : barcode_freq.data )
        {
            long num = pair.second ;
            const std::string & barcode_str = pair.first ;
            if( num > max_pair_threshold )
                barcode_list.data.AssignTag(barcode_str,0);
        }
    }

    // print barcode informations.
    void PrintBarcodeFreq()
    {
        auto bfreq = BGIQD::FILES::FileWriterFactory
            ::GenerateWriterFromFileName(fNames.barcodeFreq());
        if( bfreq == NULL )
            FATAL(" failed to open xxx.barcodeFreq to write !!! ");

        (*bfreq)<<barcode_freq.ToString();
        delete bfreq ;
    }

}config;

/**********************************************************
 *
 *  The main function
 *
 * ********************************************************/

int main(int argc , char **argv )
{
    // parse parameters ...
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , read1 , "read 1 for parse ");
        DEFINE_ARG_REQUIRED(std::string , prefix , "prefix . Output xxx.barcodeList xxx.barcodeFreq xxx.readNameList");
        DEFINE_ARG_OPTIONAL(int , min_pair_threshold , " the min_pair_threshold " , "10");
        DEFINE_ARG_OPTIONAL(int , max_pair_threshold , " the max_pair_threshold " , "500");
    END_PARSE_ARGS

    config.max_pair_threshold = max_pair_threshold.to_int() ;
    config.min_pair_threshold = min_pair_threshold.to_int() ;

    config.Init(read1.to_string() , prefix.to_string());
    BGIQD::LOG::timer timer(config.loger,"ParseReadName");

    // load read1.fasta and parse all read names.
    config.ParseRead1() ;
    // mask low quality barcode as zero-barcode.
    config.MaskTooLowBarcode() ;
    // mask too density barcode as zero-barcode because they may caused by barcode collision.
    config.MaskTooHighBarcode() ;

    // print barcode informations.
    if( config.barcode_list.preload )
    {
        config.PrintBarcodeList();
        config.PrintBarcodeFreq();
    }
}
