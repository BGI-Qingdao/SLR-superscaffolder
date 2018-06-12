#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include <map>

struct AppConf
{
    const std::string fake_read_name ;

    const std::string random_seq;

    const std::string fake_q;

    long pair_num ;
    AppConf() : 
        fake_read_name("@ST-E0:0:SIMULATE:8:0:0:")
        ,random_seq("ATCGAG")
        , fake_q(16+6,'F')
        , pair_num (0)
    {}

    std::map<int , std::string > barcodeMap ;

    void LoadBarcodeMap(std::istream & in )
    {
        auto split_2_map = [this ](const std::string & line)
        {
            auto item = BGIQD::STRING::split(line,"\t");
            if( item.size() <2 )
            {
                WARN(" invalid data from barcodeMap : "+line);
                return ;
            }
            barcodeMap[std::stoi(item[0])]=item[1];
        };
        BGIQD::FILES::FileReaderFactory::EachLine(in ,split_2_map );
    }
    void FilterAllReads(std::istream & r1 , std::istream & r2 , std::ostream & out)
    {
        std::string line1 , line2 ;
        int N= 1;
        while(!std::getline(r1,line1).eof())
        {
            auto item  = BGIQD::STRING::split(line1);
            auto itr = barcodeMap.find(std::stoi( item[1] )) ; 
            if( itr == barcodeMap.end() )
            {
                std::getline(r1,line1);
                std::getline(r1,line1);
                std::getline(r1,line1);

                std::getline(r2,line1);
                std::getline(r2,line1);
                std::getline(r2,line1);
                std::getline(r2,line1);
                continue ;
            }
            std::string barcode = itr->second ;
            out<<fake_read_name<<N<<" 1:N:0"<<std::endl;
            std::getline(r1,line1);
            out<<barcode<<random_seq<<line1<<std::endl;
            std::getline(r1,line1);
            out<<line1<<std::endl;
            std::getline(r1,line1);
            BGIQD::STRING::replace_all(line1,"!","#");
            out<<fake_q<<line1<<std::endl;

            std::getline(r2,line1);
            out<<fake_read_name<<N<<" 3:N:0"<<std::endl;
            std::getline(r2,line1);
            out<<line1<<std::endl;
            std::getline(r2,line1);
            out<<line1<<std::endl;
            std::getline(r2,line1);
            BGIQD::STRING::replace_all(line1,"!","#");
            out<<line1<<std::endl;

            pair_num++;
        }
    }

    void PrintRead2(std::ostream & out ,const std::string & indice )
    {
        for( long i = 0 ; i < pair_num ; i++ )
        {
            out<<fake_read_name<<i<<" 2:N:0"<<std::endl;
            out<<indice<<"\n+\nAAFFFKKK\n";
        }
    }
} config;



int main(int argc , char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string, map_file,"barcode map file");
    DEFINE_ARG_REQUIRED(std::string, read1 , "read1 file");
    DEFINE_ARG_REQUIRED(std::string, read2 , "read2 file");
    DEFINE_ARG_OPTIONAL(std::string, indices, "sample index","TTCACGCG");
    DEFINE_ARG_OPTIONAL(std::string, lane, "sample lane","001");
    END_PARSE_ARGS

    auto in_m = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(map_file.to_string());
    config.LoadBarcodeMap(*in_m);
    delete in_m ;


    auto in_p1 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(read1.to_string());
    auto in_p2 = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(read2.to_string());
    std::string r13f = std::string("read-RA_si_")+indices.to_string()+"_lane-"+lane.to_string()+"-chunk-001.fastq.gz";
    std::string r2f =  std::string("read-I1_si_")+indices.to_string()+"_lane-"+lane.to_string()+"-chunk-001.fastq.gz";
    auto o1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(r13f);
    config.FilterAllReads(*in_p1,*in_p2,*o1);
    delete o1 ;
    delete in_p1;
    delete in_p2;

    auto o2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(r13f);
    config.PrintRead2(*o2,indices.to_string());
    delete o2 ;

    return 0;
}
