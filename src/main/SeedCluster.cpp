/**********************************************************
 * 
 * @Brief  :
 *     Cluster candidate contigs to close a gap 
 *     based on their co-barcoding relationship.
 *
 * ********************************************************/
#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/string/stringtools.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"

#include "utils/collection/collection.h"

#include "utils/misc/fileName.h"

#include "stLFR/ContigBinBarcode.h"
#include "stLFR/ContigOrder.h"

#include <set>

//
// Struct to wrap all global variables and functions
//
struct AppConfig
{
    BGIQD::LOG::logger loger;
    BGIQD::MISC::FileNames fName;
    BGIQD::MISC::Freq<int> freqs;
    float min;
    int loop_num;
    typedef BGIQD::stLFR::OrderItem<int> GapInfo ;
    typedef BGIQD::Collection::Collection<unsigned int>  Cols;
    std::map<unsigned int , Cols> relations;
    std::set<unsigned int> trunk_seeds;
    std::set<unsigned int> PE_seeds;
    std::map<int,std::vector<GapInfo>>  infos ;
    int strategy;

    void Init( const std::string & prefix, float m , int s)
    {
        fName.Init(prefix);
        loger.Init("SeedCluster");
        min = m;
        strategy = s;
    }


    void LoadTrunk()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        BGIQD::stLFR::Load_MST_Trunk_Linear(*in,infos);
        delete in ;
        for(auto & pair: infos)
        {
            for( auto & info : pair.second )
            {
                trunk_seeds.insert(info.prev);
                trunk_seeds.insert(info.next);
            }
        }
        loger<<BGIQD::LOG::lstart() << "LoadTrunk done "<<BGIQD::LOG::lend() ;
    }

    void LoadBinCluster()
    {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fName.cluster("pe"));
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

    Cols LoopCluster(const Cols seeds)
    {
        Cols ret;
        for( auto x : seeds)
        {
            ret = Cols::Union(relations[x.first],ret);
        }
        return  Cols::Union(seeds,ret);
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
                else if( strategy == 2 )
                {
                    Cols both1 = relations[info.prev];
                    Cols both2 = relations[info.next];
                    for( int i = 1 ; i < loop_num ; i++ )
                    {
                        both1 = LoopCluster(both1);
                        both2 = LoopCluster(both2);
                    }
                    both = Cols::Intersection(both1 , both2);
                }
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
        DEFINE_ARG_REQUIRED(int, strategy, "0 for Intersection ; 1 for Union ; 2 for second cluster");
        DEFINE_ARG_OPTIONAL(int, loop_num, "loop num for strategy 2 ","1");
    END_PARSE_ARGS;
    BGIQD::LOG::timer t(config.loger,"SeedCluster");
    config.Init(prefix.to_string(), threshold.to_float(), strategy.to_int());
    config.loop_num = loop_num.to_int() ;
    config.LoadTrunk();
    config.LoadBinCluster();
    config.PrintCluster();
}
