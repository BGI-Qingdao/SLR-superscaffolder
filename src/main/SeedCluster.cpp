#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "algorithm/collection/collection.h"

#include "soap2/fileName.h"

#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <set>

struct AppConfig
{
    BGIQD::LOG::logger loger;
    BGIQD::SOAP2::FileNames fName;
    BGIQD::FREQ::Freq<int> freqs;
    float min;
    /*
    struct GapInfo
    {
        unsigned int prev ;
        unsigned int next ;
        std::set<unsigned int> relations;
    };*/
    typedef BGIQD::stLFR::TrunkGap<int> GapInfo ;
    typedef BGIQD::Collection::Collection<unsigned int>  Cols;
    std::map<unsigned int , Cols> relations;
    std::set<unsigned int> trunk_seeds;
    std::set<unsigned int> PE_seeds;
    std::map<int,std::vector<GapInfo>>  infos ;
    int strategy;

    void Init( const std::string & prefix, float m , int s)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("SeedCluster",BGIQD::LOG::loglevel::INFO, loger);
        min = m;
        strategy = s;
    }


    void LoadTrunk()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        /*
        std::string line ;
        unsigned int prev = -1 ;
        while(! std::getline(*in,line).eof() )
        {
            if( line[0] == '-' )
            {
                prev = -1 ;
                continue ;
            }
            unsigned int now = std::stoul(line);
            trunk_seeds.insert(now);
            if( prev != (unsigned int )-1 )
            {
                GapInfo info ;
                info.prev = prev ;
                info.next = now ;
                infos.push_back(info);
            }
            prev = now ;
        }
        */
        BGIQD::stLFR::Load_MST_Trunk_Linear(*in,infos);
        delete in ;
        loger<<BGIQD::LOG::lstart() << "LoadTrunk done "<<BGIQD::LOG::lend() ;
    }

    void LoadBinCluster()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.cluster());
        if( in == NULL )
            FATAL("failed to open xxx.cluster for read!!! ");
        auto parseline = [&](const std::string & line )
        {
            BGIQD::stLFR::ContigRelation tmp;
            tmp.InitFromString(line);
            if( trunk_seeds.find(tmp.contigId) == trunk_seeds.end() )
            {
                return ;
            }
            for( auto x : tmp.sims )
            {
                if( x.second.simularity >= min )
                {
                    relations[tmp.contigId].IncreaseElement(x.first);
                    PE_seeds.insert(x.first);
                }
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
        loger<<BGIQD::LOG::lstart() << "load contig sim graph done "<<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart() << "PE seeds num  "<<PE_seeds.size()<<BGIQD::LOG::lend() ;
    }

    void PrintCluster()
    {
        auto  out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.seeds_cluster_seeds());
        if( out == NULL )
            FATAL( " failed to open xxx.seeds_cluster_seeds to write"); 
        for(auto & pair: infos)
        {
            for( auto info : pair.second )
            {
                Cols both;
                if( strategy == 0 )
                    both = Cols::Intersection(relations[info.prev] , relations[info.next]);
                else if( strategy == 1 )
                    both = Cols::Union(relations[info.prev] , relations[info.next]);
                else
                    assert(0);
                (*out)<<info.prev<<'\t'<<info.next;
                freqs.Touch(both.keysize());
                for(const auto x : both)
                {
                    (*out)<<'\t'<<x.first;
                }
                (*out)<<'\n';
            }
        }
        delete out;
        loger<<BGIQD::LOG::lstart() << "PrintCluster done "<<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart() << "freq is  \n"<<freqs.ToString()<<BGIQD::LOG::lend() ;
    }

}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files.");
        DEFINE_ARG_REQUIRED(float , threshold, "min simularity threshold");
        DEFINE_ARG_REQUIRED(int, strategy, "0 for Intersection ; 1 for Union");
    END_PARSE_ARGS;
    BGIQD::LOG::timer t(config.loger,"SeedCluster");
    config.Init(prefix.to_string(), threshold.to_float(), strategy.to_int());
    config.LoadTrunk();
    config.LoadBinCluster();
    config.PrintCluster();
}
