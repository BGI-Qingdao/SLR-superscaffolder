#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"
#include "common/error/Error.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "soap2/fileName.h"

#include <vector>
#include <map>

struct AppConfig
{

    struct GapInfo
    {
        unsigned int prev ;
        unsigned int next ;
        unsigned int p1 ;
        unsigned int n1 ;
        std::vector<unsigned int> fill;
    };

    BGIQD::LOG::logger loger;
    BGIQD::SOAP2::FileNames fName;
    BGIQD::FREQ::Freq<int> freqs;

    void Init(const std::string & pre) 
    {
        fName.Init(pre);
        BGIQD::LOG::logfilter::singleton().get("FCRFByTF", BGIQD::LOG::loglevel::INFO,loger);
    }

    std::vector<GapInfo> infobuffer;
    std::map<int ,std::vector<int> > trunks;
    std::map<unsigned int , int> seed_pos ;
    std::vector< std::vector<unsigned int> > roadfill ;
    void LoadTrunk()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");

        std::string line ;
        unsigned int prev = -1 ;
        int id = 0 ;
        int pos = 0 ;
        while(! std::getline(*in,line).eof() )
        {
            if( line[0] == '-' )
            {
                prev = -1 ;
                id ++ ;
                continue ;
            }
            unsigned int now = std::stoul(line);
            //trunk_seeds.insert(now);
            if( prev != (unsigned int )-1 )
            {
                GapInfo info ;
                info.prev = prev ;
                info.next = now ;
                infobuffer.push_back(info);
                seed_pos[info.prev] =  pos ;
                trunks[id].push_back(pos);
            }
            prev = now ;
            pos ++ ;
        }
        delete in ;
        loger<<BGIQD::LOG::lstart() << "LoadTrunk done "<<BGIQD::LOG::lend() ;
    }

    void LoadTrunkFill()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.trunk_fill());
        if( in == NULL )
            FATAL(" failed to open xxx.trunk_fillfor read!!! ");
        auto loadfill = [this](const std::string & line )
        {
            unsigned int prev , next , p1 , n1 , f;
            std::istringstream ist(line); ;
            ist>>prev>>next>>p1>>n1;
            auto & gap = infobuffer[seed_pos[prev]];
            gap.p1 = p1 ;
            gap.n1 = n1 ;
            while(! ist.eof() )
            {
                ist>>f ;
                gap.fill.push_back(f);
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,loadfill);
        delete in ;
    }

    void GenerateRoadFill()
    {
        for( const auto & pair : trunks )
        {
            const auto & trunk = pair.second ;
            std::vector<unsigned int> a_fill ;
            unsigned int prev_n1 = -1 ;
            int seed_num = 1 ;
            for(size_t i = 0 ; i < trunk.size() ; i++)
            {
                const auto & info = infobuffer[trunk[i]];
                if( info.fill.empty() || ( prev_n1 != (unsigned int)-1 && prev_n1 != info.p1 ) )
                {// a conflict or unfilled gap
                    if( !a_fill.empty() )
                    {
                        roadfill.push_back(a_fill);
                        freqs.Touch(seed_num);
                    }
                    seed_num = 1 ;
                    a_fill.clear();
                    prev_n1  = -1 ;
                    continue;
                }
                a_fill.insert(a_fill.end() , info.fill.begin() , info.fill.end());
                prev_n1 = info.n1 ;
                seed_num++;
            }
        }
        loger<<BGIQD::LOG::lstart()<<"succ fill seed freq"<<freqs.ToString()<<BGIQD::LOG::lend();
    }
    void PrintRoadFill()
    {
        // fake as contig road fill
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.contigroadfill());
        for( const auto & fill : roadfill )
        {
            (*out)<<fill.size();
            for( auto i : fill)
                *(out)<<'\t'<<i;
            *(out)<<'\n';
        }
        delete out ;
    }

}config;


int  main(int argc, char **argv)
{
    //step0 Parse parmeters...
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, " In xxx.mintree_trunk_linear ; xxx.trunk_fill ");
    END_PARSE_ARGS;
    config.Init( prefix.to_string());
    config.LoadTrunk();
    config.LoadTrunkFill();
    config.GenerateRoadFill();
    config.PrintRoadFill();
}
