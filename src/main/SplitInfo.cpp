#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "biocommon/sam_bam/sam_parser.h"

#include "algorithm/distribution/distribution.h"

#include "stLFR/contigPEGraph.h"
#include "stLFR/CBB.h"
#include "stLFR/EasySam.h"

#include "soap2/fileName.h"

struct AppConfig
{

    BGIQD::LOG::logger loger;
    BGIQD::SOAP2::FileNames fName;

    BGIQD::EASY_SAM::PE_Baisc pe_basic;

    std::vector<BGIQD::EASY_SAM::PEInfo> pe_diffs;

    std::map<unsigned int , BGIQD::stLFR::ContigBarcodeInfo> cbs;

    std::map<int , BGIQD::stLFR::ContigOnBarcode> c2bs ;

    std::map<unsigned int ,BGIQD::stLFR::ContigIndex> seeds;

    BGIQD::FREQ::Freq<int> ISFreq;

    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("SplitInfo",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadContigIndex()
    {
        BGIQD::LOG::timer t(loger,"LoadContigIndex");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.ContigIndex()) ;
        if( in == NULL )
            FATAL( "open .ContigIndex file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            BGIQD::stLFR::ContigIndex tmp;
            tmp.InitFromString(line);
            seeds[tmp.contig] = tmp;
        }
        delete in ;
    };

    void ParseRead2Contig()
    {
        BGIQD::LOG::timer t(loger,"ParseRead2Contig");
        std::string line ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.read2contig());
        if ( in == NULL )
            FATAL(" open xxx.read2contig to read failed !!! ");

        BGIQD::EASY_SAM::EasySam prev_line;
        prev_line.contig_name = 0 ;

        while( ! std::getline( *in , line ).eof() )
        {
            BGIQD::EASY_SAM::EasySam curr_line;
            curr_line.InitFromString(line);
            if( curr_line.barcode != 0 && parse_barcode )
            {
                cbs[curr_line.contig_name].Touch(curr_line.pos_1bp,curr_line.barcode);
                c2bs[curr_line.barcode].Touch(curr_line.contig_name);
            }

            if( curr_line.pe_match && parse_pe )
            {
                if( (curr_line.is_p ||  (! curr_line.pe_match) ))
                {
                    prev_line = curr_line ;
                    continue;
                }
                if( !( prev_line.is_p && prev_line.pe_match) ) 
                {
                    prev_line = curr_line ;
                    continue;
                }
                if( prev_line.contig_name == curr_line.contig_name)
                {
                    std::vector<int> pos;
                        pos.push_back(prev_line.pos_1bp);
                        pos.push_back(curr_line.pos_1bp);
                    if( prev_line.match_reverse )
                        pos.push_back(prev_line.pos_1bp-99);
                    else
                        pos.push_back(prev_line.pos_1bp+99);
                    if( curr_line.match_reverse)
                        pos.push_back(curr_line.pos_1bp-99);
                    else
                        pos.push_back(curr_line.pos_1bp+99);

                    std::sort(pos.begin() ,pos.end());
                    int IS = pos[3] - pos [0] +1 ;
                    ISFreq.Touch(IS);
                }
                else
                {
                    BGIQD::EASY_SAM::PEInfo tmp ;
                    tmp.read1 = prev_line.read_id ;
                    tmp.read2 = curr_line.read_id ;
                    tmp.match_reverse1 = prev_line.match_reverse ;
                    tmp.match_reverse2 = curr_line.match_reverse ;
                    tmp.pos_1bp1 = prev_line.pos_1bp; 
                    tmp.pos_1bp2 = curr_line.pos_1bp;
                    tmp.contig1 = prev_line.contig_name ;
                    tmp.contig2 = curr_line.contig_name;
                    pe_diffs.push_back(tmp);
                }
            }
            prev_line = curr_line ;
        }
        for( auto & pair : cbs )
        {
            pair.second.contig_id = pair.first ;
        }
        for( auto & pair : c2bs )
        {
            pair.second.barcode_id = pair.first ;
        }
    }

    bool parse_pe;
    bool parse_barcode;

    bool CheckParameters() const 
    {
        return parse_pe || parse_barcode ;
    }

    void PrintResult()
    {
        if( parse_pe) 
        {
            auto in1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.pe_pairs());
            for( const auto & item : pe_diffs )
            {
                (*in1)<<item.ToString()<<'\n';
            }
            delete in1 ;

            auto in2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.pe_info());
            long total = 0;
            long num = 0;
            for( const auto & pair : ISFreq.data)
            {
                total += pair.first * pair.second;
                num += pair.second ;
            }
            (*in2)<<"Average insert size is :\t"<<total/num<<'\n';
            delete in2;
        }
        if(parse_barcode)
        {
            auto in1 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.BarcodeOnContig());
            for( const auto & pair : cbs)
            {
                (*in1)<<pair.second.ToString()<<'\n';
            }
            delete in1 ;
            auto in2 = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.contigOnBarcode());
            for( const auto & pair : c2bs)
            {
                (*in2)<<pair.second.ToString()<<'\n';
            }
            delete in2 ;
        }
    }

}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files. Input xxx.read2contig , Output depends on below options. Please at least choose one of below options.");
        DEFINE_ARG_OPTIONAL(bool, parse_pe,"parse pe data. will output xxx.pe_info && xxx.pe_pairs", "false");
        DEFINE_ARG_OPTIONAL(bool, parse_barcode,"parse barcode data . will output xxx.barcodeOnContig && xxx.contigOnBarcode", "false");
    END_PARSE_ARGS;

    config.parse_pe = parse_pe.to_bool() ;
    config.parse_barcode = parse_barcode.to_bool() ;
    config.Init(prefix.to_string());
    if( ! config.CheckParameters()  )
    {
        config.loger<<BGIQD::LOG::lstart()<<"Nothing to do ... Exit now."<<BGIQD::LOG::lend();
        return 0;
    }

    BGIQD::LOG::timer t(config.loger,"SplitInfo");
    config.LoadContigIndex();
    config.ParseRead2Contig();
    config.PrintResult();
    return 0 ;
}
