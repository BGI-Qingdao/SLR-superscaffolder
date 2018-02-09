#include "soap2/contigGraph.h"
#include "soap2/loadGraph.h"
#include <iostream>
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

void report(const  BGIQD::SOAP2::GlobalConfig & config)
{
    std::map<int,int> freq;
    std::cout<<"--- Paths start  ----"<<std::endl;
    for(const auto & i : config.contigs)
    {
        if( freq.find(i.size()) == freq.end() )
            freq[i.size()] = 1;
        else
            freq[i.size()] ++;
        for( auto j : i)
            std::cout<<j<<'\t';
        std::cout<<std::endl;
    }
    std::cout<<"--- Paths end ----"<<std::endl;

    std::cout<<"clusterNum "<<config.clusterNum<<std::endl;
    std::cout<<"pathNum "<<config.contigs.size()<<std::endl;
    std::cout<<"--- freq start ---- "<<std::endl;
    for( const auto & i : freq)
        std::cout<<i.first<<"\t"<<i.second<<std::endl;
    std::cout<<"--- freq end ---- "<<std::endl;
}

int main(int argc , char **argv)
{
    BGIQD::LOG::logger lger;
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , lger);
    BGIQD::LOG::timer t(lger,"SuperContig");
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix");
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    END_PARSE_ARGS
    lger<<BGIQD::LOG::lstart()<<"parse args end ... "<<BGIQD::LOG::lend();

    BGIQD::SOAP2::GlobalConfig config;
    config.K = kvalue.to_int();
    config.arc = prefix.to_string() +".Arc";
    config.updateEdge = prefix.to_string() +".updated.edge";
    config.cluster= prefix.to_string() +".cluster";

    lger<<BGIQD::LOG::lstart()<<"loadUpdateEdge start ... "<<BGIQD::LOG::lend();
    BGIQD::SOAP2::loadUpdateEdge(config);
    lger<<BGIQD::LOG::lstart()<<"loadArc start ... "<<BGIQD::LOG::lend();
    BGIQD::SOAP2::loadArc(config);
    lger<<BGIQD::LOG::lstart()<<"loadCluster start ... "<<BGIQD::LOG::lend();
    BGIQD::SOAP2::loadCluster(config);
    lger<<BGIQD::LOG::lstart()<<"buildConnection start ... "<<BGIQD::LOG::lend();
    BGIQD::SOAP2::buildConnection(config);
    lger<<BGIQD::LOG::lstart()<<"report start ... "<<BGIQD::LOG::lend();
    report(config);
    return 0;
}
