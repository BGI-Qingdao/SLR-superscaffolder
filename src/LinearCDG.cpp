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
    std::vector<BGIQD::SOAP2::ContigRoad> contigs ; 
    void Init()
    {
        BGIQD::LOG::logfilter::singleton().get("LinearCDG",BGIQD::LOG::loglevel::INFO , loger );
    }

}config;

int main(int argc , char ** argv)
{

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(bool, solve, 's',true,"solve multi ? default not");
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
                index += 2 ;
                linear_len.push_back( curr.GetValidTo().length );
                linear_sim.push_back( curr.GetValidTo().sim);
                linear_len.push_back( curr.GetValidFrom().length );
                linear_sim.push_back( curr.GetValidFrom().sim);
            }
            /*
            else if ( curr.IsTipTo()  && curr.to_size == 1 )
            {
                index ++ ;
                linear_len.push_back( curr.GetValidTo().length );
                linear_sim.push_back( curr.GetValidTo().sim);
            }*/
        }
    }
    config.loger<<BGIQD::LOG::lstart()<<" total linear node "<<index<<BGIQD::LOG::lend();
    std::sort(linear_len.begin() , linear_len.end());
    std::sort(linear_sim.rbegin() , linear_sim.rend());
    for(int i = 0 ; i<index ; i++ )
    {
        config.loger<<BGIQD::LOG::lstart()<<linear_len[i]<<"\t"<<linear_sim[i]<<BGIQD::LOG::lend(); 
    }
    int len_threshold = linear_len[ index * 0.95 ];
    float sim_threshold = linear_sim[ index * 0.95 ];
    config.loger<<BGIQD::LOG::lstart()<<" used linear len threshold "<<len_threshold<<BGIQD::LOG::lend();
    config.loger<<BGIQD::LOG::lstart()<<" used linear sim threshold "<<sim_threshold<<BGIQD::LOG::lend();
    // solve multi

    for( auto & m : config.edges )
    {
        auto & curr = m.second;
        if( curr.IsLinear() 
                || (curr.IsTipTo() && curr.to_size == 1 ) 
                || (curr.IsTipFrom() && curr.from_size == 1 ) )
        {
            continue ;
        }
        if( curr.from_size > 1 )
        {
            for( auto & i : curr.from )
            {
                if( i.second.length > len_threshold || i.second.sim > sim_threshold )
                    i.second.SetJump();
            }
        }
        if( curr.to_size > 1 )
        {
            for( auto & i : curr.to)
            {
                if( i.second.length > len_threshold || i.second.sim > sim_threshold )
                    i.second.SetJump();
            }
        }
    }
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
    BGIQD::FREQ::Freq<int> len_freq;

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

    for(auto & m : config.edges)
    {
        auto & curr = m.second;

        if( curr.IsMarked() 
                ||curr.IsLinear() || curr.IsCircle()
          )
            continue ;
        if( curr.IsSingle() )
        {
            curr.Mark();
        }
        else
        {
            BGIQD::SOAP2::ContigRoad path;
            auto extractPath = [&curr,&path]( unsigned int to_s , bool to_order, bool search_order)
            {
                unsigned int to = to_s ;
                path.contig.clear();
                path.real_contig.clear();
                if( search_order )
                {
                    path.downstream = true ;
                    path.real_contig.push_back(curr.edge_id);
                }
                else
                {
                    path.downstream = false ;
                    path.real_contig.push_back(curr.bal_id);
                }
                //detect headin
                if( search_order && curr.to_size == 1 )
                {
                    path.headin = true;
                }
                else if ( ! search_order && curr.from_size == 1 )
                {
                    path.headin = true;
                }
                else
                {
                    path.headin = false ;
                }
                path.contig.push_back(curr.edge_id) ;
                path.contig.push_back(to) ;

                if( search_order ^ to_order )
                {
                    path.real_contig.push_back(config.edges.at(to).bal_id);
                }
                else
                {
                    path.real_contig.push_back(config.edges.at(to).edge_id);
                }
                bool order = search_order;//
                bool torder = to_order; // 
                if( config.edges.at(to).IsMarked())
                    return ;
                while( (! config.edges.at(to).IsCircle())&& config.edges.at(to).IsLinear() && !config.edges.at(to).IsMarked() )
                {
                    auto & curr_to = config.edges.at(to);
                    curr_to.Mark();
                    //downstream
                    if( order )
                    {
                        //
                        //A1->B1
                        if( torder )
                        {
                            const auto & conn  =  curr_to.GetValidTo() ;
                            to = conn.to;
                            torder = conn.IsPositive();
                            if( torder )
                            {
                                path.real_contig.push_back(config.edges.at(to).edge_id );
                            }
                            else
                            {
                                path.real_contig.push_back(config.edges.at(to).bal_id);
                            }
                        }
                        //A1->B2
                        else
                        {
                            const auto & conn  =  curr_to.GetValidFrom() ;
                            to = conn.to;
                            torder = conn.IsPositive();
                            if( torder )
                            {
                                path.real_contig.push_back(config.edges.at(to).bal_id);
                            }
                            else
                            {
                                path.real_contig.push_back(config.edges.at(to).edge_id);
                            }
                            order = false;
                        }
                    }
                    //upstream
                    else
                    {
                        //A1<-B1
                        //A2->B2
                        if( torder )
                        {

                            const auto & conn  =  curr_to.GetValidFrom() ;
                            to = conn.to;
                            torder = conn.IsPositive();
                            if( torder )
                            {
                                path.real_contig.push_back(config.edges.at(to).bal_id);
                            }
                            else
                            {
                                path.real_contig.push_back(config.edges.at(to).edge_id);
                            }
                            order = false;
                        }
                        //A1<-B2
                        //A2->B1
                        else
                        {
                            const auto & conn  =  curr_to.GetValidTo() ;
                            to = conn.to;
                            torder = conn.IsPositive();
                            if( torder )
                            {
                                path.real_contig.push_back(config.edges.at(to).edge_id );
                            }
                            else
                            {
                                path.real_contig.push_back(config.edges.at(to).bal_id);
                            }
                        }
                    }
                    path.contig.push_back(to) ;
                }
                if( config.edges.at(to).IsMarked() )
                    return ;
                // line->next_k->
                // <-next_k<-line
                path.tailin = false;

                auto & curr_to = config.edges.at(to);
                if( ! curr_to.IsCircle() )
                {

                    if( ( torder && order ) || ( !torder && !order ) )
                    {
                        if( curr_to.from_size == 1 )
                        {
                            path.tailin = true;
                        }
                    }
                    else
                    {
                        if( curr_to.to_size == 1 )
                        {
                            path.tailin = true;
                        }
                    }
                }
                path.length = path.contig.size();
                if(! path.headin )
                    path.length --;
                if(!path.tailin )
                    path.length -- ;
                {
                    config.contigs.push_back(path);
                }
            };
            for( auto next : curr.to )
            {
                if(! next.second.IsValid() )
                    continue;
                extractPath(next.second.to,next.second.IsPositive(),true);
            }
            for(auto next : curr.from)
            {
                if( ! next.second.IsValid() )
                    continue;
                extractPath(next.second.to,next.second.IsPositive(),false);
            }
            curr.Mark();
        }
    }

    for(const auto & i : config.contigs)
    {
        if( i.length < 2 )
            continue;
        std::cout<<i.length<<'\t';
        int start = 1 ;
        if( i.headin )
        {
            start --;
        }
        len_freq.Touch(i.length);
        for( int j =0 ; j<i.length ; j++ )
        {
            std::cout<<i.real_contig[start+j]<<"\t";
        }
        std::cout<<std::endl;
    }
    config.loger<<BGIQD::LOG::lstart()<<"len freq"<< '\n'<<len_freq.ToString() << BGIQD::LOG::lend();
    return 0 ;
}
