#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"

#include "soap2/fileName.h"
#include "stLFR/EasySam.h"

#include "algorithm/incr_array/incr_array.h"
#include "algorithm/multi_key_hash/MultiKeyHash.h"

#include <iostream>
#include <set>
#include <string>

struct AppConfig
{

    typedef BGIQD::MultiKeyMap::
        BiKeyHash<unsigned int , int > ConnMap;

    typedef BGIQD::EASY_SAM::EasySam_V1 SAMInfo;

    typedef BGIQD::INCRARRAY::IncrArray
        <SAMInfo> SAMLoadBuffer;


    struct ReadPair
    {
        std::set<const SAMInfo *> r1;
        std::set<const SAMInfo *> r2;
    };

    std::map< long , ReadPair > read_pair_buffer;


    ConnMap FRConn;

    ConnMap FFConn ;

    SAMLoadBuffer load_buffer ;

    BGIQD::LOG::logger loger ;

    BGIQD::SOAP2::FileNames fName ;

    void Init( const std::string & prefix)
    {
        BGIQD::LOG::logfilter::singleton().get("MergePEInfo",BGIQD::LOG::DEBUG,loger);
        fName.Init(prefix);
        load_buffer.Init(10000);
    }

    void LoadEasySamLine(const std::string & line)
    {
        BGIQD::EASY_SAM::EasySam_V1 tmp ;
        tmp.InitFromString(line);
        load_buffer.push_back(tmp);
    }

    void LoadRead2Contig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(
                fName.read2contig_v1()
                );
        if( in == NULL )
            FATAL("failed to open xxx.read2contig_v1 for read !!!");
        BGIQD::FILES::FileReaderFactory::EachLine(*in
                ,std::bind(
                    &AppConfig::LoadEasySamLine,this,std::placeholders::_1)
                );
        delete in ;
    }

    void LoadContig2Read()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(
                fName.contig2read_v1()
                );
        if( in == NULL )
            FATAL("failed to open xxx.contig2read_v1 for read !!!");
        BGIQD::FILES::FileReaderFactory::EachLine(*in
                ,std::bind(
                    &AppConfig::LoadEasySamLine,this,std::placeholders::_1)
                );
        delete in ;
    }

    void BuildConns()
    {
        for(const SAMInfo & item : load_buffer )
        {
            auto r = item.read_id ;
            if( item.read_index == 1 )
                (read_pair_buffer[r].r1).insert(&item);
            else
                (read_pair_buffer[r].r2).insert(&item);
        }

        for( const auto & pair : read_pair_buffer )
        {
            for( const auto & r1 : pair.second.r1 )
            {
                const auto & info1 = *r1 ;
                for( const auto & r2 : pair.second.r2 )
                {
                    const auto & info2 = *r2 ;
                    if( info1.match_reverse ^ info2.match_reverse )
                    {
                        if( FFConn.Contain( info1.contig_name , info2.contig_name ) )
                        {
                            FFConn.Set( info1.contig_name , info2.contig_name 
                                    , FFConn.At( info1.contig_name , info2.contig_name ) + 1);
                        }
                        else
                        {
                            FFConn.Set( info1.contig_name , info2.contig_name 
                                    , 1 );

                        }
                    }
                    else
                    {
                        if( FRConn.Contain( info1.contig_name , info2.contig_name ) )
                        {
                            FRConn.Set( info1.contig_name , info2.contig_name 
                                    , FRConn.At( info1.contig_name , info2.contig_name ) + 1);
                        }
                        else
                        {
                            FRConn.Set( info1.contig_name , info2.contig_name 
                                    , 1 );

                        }
                    }
                }
            }
        }
    }

    void PrintConns()
    {
        auto out = BGIQD::FILES::FileWriterFactory::
            GenerateWriterFromFileName(fName.contig_pe_conns());
        if( out == NULL )
            FATAL("failed to open xxx.contig_pe_conns for write !!!");

        for( const auto & item : FFConn.data)
        {
            (*out)<<item.first.first<<'\t'
                <<item.first.second<<'\t'
                <<item.second<<'\t'<<'+'<<'\n';
        }
        for( const auto & item : FRConn.data)
        {
            (*out)<<item.first.first<<'\t'
                <<item.first.second<<'\t'
                <<item.second<<'\t'<<'-'<<'\n';
        }
        delete out;
    }
} config;

int main(int argc , char ** argv)
{
    // parse args
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string,prefix, "prefix. \n\
                                                Input\n\
                                                    xxx.read2contig_v1\n\
                                                    xxx.contig2read_v1\n\
                                                Output xxx.pe_pair_v1");
    END_PARSE_ARGS

    config.Init(prefix.to_string());

    BGIQD::LOG::timer t(config.loger,"ParseRead2Contig");
    config.LoadRead2Contig();
    config.LoadContig2Read();
    config.BuildConns();
    config.PrintConns();
    return 0;
}
