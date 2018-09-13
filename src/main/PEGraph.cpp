#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "stLFR/contigPEGraph.h"
#include "stLFR/EasySam.h"
#include "soap2/fileName.h"

struct AppConfig
{
    BGIQD::LOG::logger loger;
    BGIQD::stLFR::ContigPEGraph pe_graph;
    BGIQD::SOAP2::FileNames fName;

    std::map<unsigned int , int > contigLen_cache ;
    std::map<unsigned int , std::map<unsigned int , std::vector<int> > > pe_cache ;

    int insert_size ;
    int max_is ;

    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            auto items = BGIQD::STRING::split(line,"\t");
            contigLen_cache[std::stoul(items[0])] = std::stoul(items[1]);
            contigLen_cache[std::stoul(items[0])+1] = std::stoul(items[1]);
            pe_graph.AddNode(std::stoul(items[0]) , std::stoul(items[1]));
        }
        delete in ;
    };



    void LoadPEInfo()
    {
        std::string line ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.pe_info());
        if ( in == NULL )
            FATAL(" open xxx.pe_info to read failed !!! ");
        std::getline(*in , line);
        delete in ;
        sscanf(line.c_str() , "Aveage insert size is :\t%d",&insert_size);
    }

    void LoadPEPair()
    {
        BGIQD::LOG::timer t(loger,"LoadPECahce");
        BGIQD::FREQ::Freq<int> ISFreq;
        std::string line ;

        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.pe_pairs());
        if ( in == NULL )
            FATAL(" open xxx.pe_pair to read failed !!! ");

        while( ! std::getline( *in , line ).eof() )
        {
            BGIQD::EASY_SAM::PEInfo tmp ;
            tmp.InitFromString(line);

            unsigned int R1Contig = tmp.contig1 , R2Contig = tmp.contig2;
            unsigned int R1Contig1 = tmp.contig1+1 , R2Contig1 = tmp.contig2+1;
            int R1CLeft = 0, R1CRight =0 , R2CLeft =0, R2CRight =0;

            //Zero at R1's left most 1 bp
            if(  tmp.match_reverse1 )
            {
                R1Contig ++ ;
                R1Contig1 -- ;
                R1CRight =  ( tmp.pos_1bp1 - 1 ) ;
                R1CLeft = -( contigLen_cache[tmp.contig1] - tmp.pos_1bp1 ) ;
            }
            else
            {
                R1CLeft = - tmp.pos_1bp1 +1 ;
                R1CRight = contigLen_cache[tmp.contig1] - tmp.pos_1bp1 ;
            }
            if( ! tmp.match_reverse2 )
            {
                R2Contig ++ ;
                R2Contig1 -- ;
                R2CLeft = (insert_size -1) - ( contigLen_cache[tmp.contig2] - tmp.pos_1bp2 ) ;
                R2CRight = (insert_size -1) + tmp.pos_1bp2 -1 ;
            }
            else
            {
                R2CLeft = (insert_size -1) - (tmp.pos_1bp2 -1 ) ;
                R2CRight = (insert_size -1) + contigLen_cache[tmp.contig2] - tmp.pos_1bp2;
            }

            if( R1CLeft < R2CLeft && R1CRight < R2CRight ) 
            { // R1C -> R2C
                int  gap = R2CLeft - R1CRight ;
                if( insert_size - gap > max_is ) 
                    continue ;
                pe_cache[R1Contig][R2Contig].push_back(gap);
                pe_cache[R2Contig1][R1Contig1].push_back(gap);
            }
            else if ( R1CLeft > R2CLeft && R1CRight > R2CRight )
            { // R2C -> R1C
                int  gap = R1CLeft - R2CRight ;
                if( insert_size - gap > max_is ) 
                    continue ;
                pe_cache[R2Contig][R1Contig].push_back(gap);
                pe_cache[R1Contig1][R2Contig1].push_back(gap);
            }
            else
            {
                //TODO
            }
        }
        delete in ;
    }

    void BuildPEGraph()
    {
        for( const auto & pair : pe_cache )
        {
            unsigned int PcontigId = pair.first ;
            for( const auto & pair1 : pair.second )
            {
                unsigned int EcontigId = pair1.first ;
                int total_len = 0 ;
                for( auto left : pair1.second )
                {
                    total_len += left ;
                }
                int average_len = total_len / (int) pair1.second.size() ;
                pe_graph.AddEdge(PcontigId ,EcontigId , average_len , pair1.second.size());
            }
        }
    }

    void PrintContigPEGraph()
    {
        BGIQD::LOG::timer t(loger,"PrintContigPEGraph");
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.pe_graph());
        if( out == NULL )
            FATAL( " failed to open xxx.pe_graph for write ");
        pe_graph.PrintAsDOT(*out);
        delete out;
    }

}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of read name \n\
                                                    Input \n\
                                                        xxx.seeds\n\
                                                        xxx.pe_info\n\
                                                        xxx.pe_pairs\n\
                                                    Output \n\
                                                        xxx.pe_graph");
    DEFINE_ARG_OPTIONAL(int, max_is ,"max valid insert_size","1000");
    END_PARSE_ARGS;

    config.max_is= max_is.to_int();
    config.Init(prefix.to_string());
    BGIQD::LOG::timer t(config.loger,"PEGraph");

    config.LoadSeeds();
    config.LoadPEInfo();
    config.LoadPEPair();
    config.BuildPEGraph();
    config.PrintContigPEGraph();
    return 0 ;
}
