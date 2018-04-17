#include <cassert>
#include <iostream>
#include <map>
#include <tuple>
#include <algorithm>
#include <vector>
#include <set>
#include "common/string/stringtools.h"
#include "common/args/argsparser.h"
#include "common/freq/freq.h"
#include "soap2/contigGraph.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

struct  AppGlobalData
{
    std::map<unsigned int , BGIQD::SOAP2::KeyEdge> edges;
    BGIQD::LOG::logger loger;

    void Init()
    {
        BGIQD::LOG::logfilter::singleton().get("LinearCDG",BGIQD::LOG::loglevel::INFO , loger );
    }

}config;

int main(int argc , char ** argv)
{

    START_PARSE_ARGS
        DEFINE_ARG_DETAIL(int , type, 'o',false,"1 for len , 2 for simularity \
                ,3 for split multi and base \
                ,4 for solve multi");
    DEFINE_ARG_DETAIL(int , len_max, 'a',false,"max_len for first close");
    DEFINE_ARG_DETAIL(int , len_min, 'i',false,"min len for second close");
    DEFINE_ARG_DETAIL(float, sim, 's',false,"min simularity");
    END_PARSE_ARGS

    config.Init();
        // Load graph
        std::string line ;
    int id  = 0 ;
    while( ! std::getline( std::cin , line).eof() )
    {
        auto items1 = BGIQD::STRING::split( line , "\t") ;
        assert( items1.size() >= 3 );
        unsigned int key = std::stoul(items1[0]);

        if( config.edges.find( key ) == config.edges.end() )
        {
            config.edges[key].Init( ++id , key );
        }

        if( items1[1] == "+" )
        {
            for(size_t i = 2 ; i < items1.size() ; i++ )
            {
                config.edges.at(key).InitTo(items1[i]);
            }
        }
        else
        {
            for(size_t i = 2 ; i < items1.size() ; i++ )
            {
                config.edges.at(key).InitFrom(items1[i]);
            }
        }
    }

    // report before action

    // anlysis already linear part

    for( auto & m : config.edges )
    {
        m.second.SetType();
    }

    std::vector<int> linear_len;

    std::vector<float> linear_sim;

    int index = 0 ;
    for( auto & m : config.edges )
    {
        auto & curr = m.second;
        if( ! curr.IsCircle() )
        {
            if( curr.IsLinear() )
            {
                index ++ ;
                linear_len.push_back( curr.GetValidTo().length );
                linear_sim.push_back( curr.GetValidTo().sim);
                linear_len.push_back( curr.GetValidFrom().length );
                linear_sim.push_back( curr.GetValidFrom().sim);
            }
        }
    }
    config.loger<<BGIQD::LOG::lstart()<<" total linear node "<<index<<BGIQD::LOG::lend();
    std::sort(linear_len.begin() , linear_len.end());
    std::sort(linear_sim.rbegin() , linear_sim.rend());

    int len_threshold = linear_len[ index * 0.8 ];
    int sim_threshold = linear_sim[ index * 0.9 ];
    config.loger<<BGIQD::LOG::lstart()<<" used linear len threshold "<<len_threshold<<BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<" used linear sim threshold "<<sim_threshold<<BGIQD::LOG::lend();
    // report after action

    for( auto & m : config.edges)
    {
        m.second.SetType();
    }
    BGIQD::FREQ::Freq<std::string> freq;
    BGIQD::FREQ::Freq<std::string> basefreq;
    BGIQD::FREQ::Freq<int> from;
    BGIQD::FREQ::Freq<int> to;
    BGIQD::FREQ::Freq<int> total;
    BGIQD::FREQ::Freq<int> del;
    BGIQD::FREQ::Freq<int> base_from;
    BGIQD::FREQ::Freq<int> base_to;
    BGIQD::FREQ::Freq<int> in_circle;

    for( const auto & m : config.edges)
    {

        auto & curr = m.second;
        if( curr.from.size() == 0 && curr.to.size() == 0 )
            basefreq.Touch("Single");
        else if( curr.from.size() == 0 && curr.to.size() == 1 )
            basefreq.Touch("Tipto 1");
        else if( curr.from.size() == 0 && curr.to.size() > 1 )
            basefreq.Touch("Tipto >1");
        else if( curr.from.size() == 1 && curr.to.size() ==0 )
            basefreq.Touch("Tipfrom 1");
        else if( curr.from.size() > 1 && curr.to.size() ==0 )
            basefreq.Touch("Tipfrom >1");
        else if( curr.from.size() == 1 && curr.to.size() == 1 )
            basefreq.Touch("Linear");
        else
            basefreq.Touch("Multi");


        if( curr.IsSingle() )
            freq.Touch("Single");
        else if ( curr.IsLinear() )
            freq.Touch("Linear");
        else if ( curr.IsTipTo() && curr.to_size ==1 )
            freq.Touch("Tipto 1");
        else if ( curr.IsTipTo() && curr.to_size > 1 )
            freq.Touch("Tipto >1");
        else if ( curr.IsTipFrom() && curr.from_size ==1 )
            freq.Touch("Tipfrom 1");
        else if ( curr.IsTipFrom() && curr.from_size > 1 )
            freq.Touch("Tipfrom >1");
        else
            freq.Touch("Multi");
        from.Touch(curr.from_size);
        to.Touch(curr.to_size);
        total.Touch(curr.total_size);
        del.Touch(curr.jump_conn);
        base_from.Touch( curr.from.size());
        base_to.Touch( curr.to.size());
        in_circle.Touch(curr.IsCircle());
    }

    config.loger<<BGIQD::LOG::lstart()<<"base key type freq"<<'\n'<<basefreq.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"key type freq"<<'\n'<< freq.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"from freq"<< '\n'<<from.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"to freq"<< '\n'<<to.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"total freq"<< '\n'<<total.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"delete freq"<< '\n'<<del.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"basefrom freq"<< '\n'<<base_from.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"baseto freq"<< '\n'<<base_to.ToString() << BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<"incircle freq"<< '\n'<<in_circle.ToString() << BGIQD::LOG::lend();

    // print linear collection

    return 0 ;
}
