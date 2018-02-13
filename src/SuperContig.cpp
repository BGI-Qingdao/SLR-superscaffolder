#include "soap2/contigGraph.h"
#include "soap2/loadGraph.h"
#include <iostream>
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/freq/freq.h"
#include <atomic>

BGIQD::LOG::logger lger;
std::atomic<int> index ;


void findConnection(BGIQD::SOAP2::GlobalConfig & config
        , unsigned int edge_id
        , bool order
        )
{
    BGIQD::SOAP2::Edge & root = config.edge_array[edge_id];
    unsigned int i = edge_id;
    std::stack<BGIQD::SOAP2::Edge> stack;
    std::map<unsigned int , BGIQD::SOAP2::Edge > history;
    std::map<unsigned int , std::vector<std::stack<BGIQD::SOAP2::Edge> > > paths;
    std::map<unsigned int , std::vector<std::stack<BGIQD::SOAP2::Edge> > > mids;
    if ( order )
        config.edge_array[i].DepthSearch( config.edge_array , stack,
            history, paths , mids ,config.edge_array[i].length , config.connections.at(root.bal_id) );
    else 
        config.edge_array[i].DepthSearch( config.edge_array , stack,
            history, paths , mids ,config.edge_array[i].length , config.connections.at(root.id) );

        index ++ ;
        for(auto & j : paths)
        {
            /*std::cerr<<j.first<<'\t'<<j.second[0].size()<<'\t';
            while(j.second[0].size() > 0 )
            {
                std::cerr<<j.second[0].top().id<<"\t";
                j.second[0].pop();
            }
            std::cerr<<std::endl;*/
            unsigned int to_id = j.first ;
            unsigned int to_id_in_path = j.second[0].top().id ;
            //
            //  A1->B1
            //  A2<-B2
            //
            if( !order && to_id == to_id_in_path)
            {
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[i]]);
                    config.key_array[config.key_map[i]].to[j.first]= true;
                }
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                    config.key_array[config.key_map[j.first]].from[(i)] = true ;
                }
            }
            //
            // A1->B2
            // A2<-B1
            //
            else if( !order && to_id != to_id_in_path )
            {
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[i]]);
                    config.key_array[config.key_map[i]].to[j.first]= false;
                }
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                    config.key_array[config.key_map[j.first]].to[(i)] = false;
                }
            }
            //
            // A2->B1
            // A1<-B2
            //
            else if ( order && to_id == to_id_in_path )
            {
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[root.bal_id]]);
                    config.key_array[config.key_map[root.bal_id]].from[j.first]= false;
                }
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                    config.key_array[config.key_map[j.first]].from[(root.bal_id)] = false;
                }
            }
            //
            // A2->B2
            // A1<-B1
            //
            else
            {
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[root.bal_id]]);
                    config.key_array[config.key_map[root.bal_id]].from[j.first]= true;
                }
                {
                    std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                    config.key_array[config.key_map[j.first]].to[(root.bal_id)] = true;
                }
            }
        }
        if( index %100 == 0 )
        {
            std::lock_guard<std::mutex> lm(config.contig_mutex);
            lger<<BGIQD::LOG::lstart()<<"process "<<index<<" ..."<<BGIQD::LOG::lend();
        }
}

void linearConnection(BGIQD::SOAP2::GlobalConfig &config , unsigned int key_id)// , BGIQD::MultiThread::MultiThread & queue)
{

    unsigned int i = key_id;
    BGIQD::SOAP2::KeyEdge & curr =  config.key_array[i];
    if( curr.IsMarked() 
            ||curr.IsLinear() 
      )
        return ;
    index ++ ;
    if( curr.IsSingle() )
    {
        std::vector<unsigned int > path;
        path.push_back(curr.edge_id) ;
        {
            std::lock_guard<std::mutex> lm(config.contig_mutex);
            config.contigs.push_back(path);
        }
        curr.Mark();
    }
    else
    {

        std::vector<unsigned int > path;
        auto extractPath = [&config,&curr,&path]( unsigned int to , bool to_order, bool search_order)
        {
            path.clear();
            path.push_back(curr.edge_id) ;
            path.push_back(to) ;
            unsigned int next_k = config.key_map[to];
            bool order = search_order;
            bool torder = to_order;
            if( config.key_array[next_k].IsMarked())
                return ;
            while( config.key_array[next_k].IsLinear() )
            {
                config.key_array[next_k].Mark();
                unsigned int next_i;
                //downstream
                if( order )
                {
                    //
                    //A1->B1
                    if( torder )
                    {
                        auto itr = config.key_array[next_k].to.begin() ;
                        next_i = itr->first;
                        torder = itr->second;
                    }
                    //A1->B2
                    else
                    {
                        auto itr = config.key_array[next_k].from.begin() ;
                        next_i = itr->first;
                        torder = itr->second;
                        order = false;
                    }
                }
                //upstream
                else
                {
                    //A1<-B1
                    if( torder )
                    {
                        auto itr = config.key_array[next_k].from.begin() ;
                        next_i = itr->first;
                        torder = itr->second;
                    }
                    //A1<-B2
                    else
                    {
                        auto itr = config.key_array[next_k].to.begin() ;
                        next_i = itr->first;
                        torder = itr->second;
                        order = true;
                    }
                }
                next_k = config.key_map[next_i];
                path.push_back(next_i) ;
            }

            {
                std::lock_guard<std::mutex> lm(config.contig_mutex);
                config.contigs.push_back(path);
            }

        };
        for(auto next : curr.to )
        {
            extractPath(next.first,next.second,true);
        }
        for(auto next : curr.from)
        {
            extractPath(next.first,next.second,false);
        }
        curr.Mark();
    }
    if( index %100 == 0 )
    {
        std::lock_guard<std::mutex> lm(config.contig_mutex);
        lger<<BGIQD::LOG::lstart()<<"process "<<index<<" ..."<<BGIQD::LOG::lend();
    }
}

void report(const  BGIQD::SOAP2::GlobalConfig & config)
{
    BGIQD::FREQ::Freq<unsigned int > freq;
    std::cout<<"--- Paths start  ----"<<std::endl;
    for(const auto & i : config.contigs)
    {
        freq.Touch(i.size());
        for( auto j : i)
            std::cout<<j<<'\t';
        std::cout<<std::endl;
    }
    std::cout<<"--- Paths end ----"<<std::endl;

    std::cout<<"clusterNum "<<config.clusterNum<<std::endl;
    std::cout<<"pathNum "<<config.contigs.size()<<std::endl;
    std::cout<<"--- freq start ---- "<<std::endl;
        std::cout<<freq.ToString();
    std::cout<<"--- freq end ---- "<<std::endl;
}

int main(int argc , char **argv)
{
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , lger);
    BGIQD::LOG::timer t(lger,"SuperContig");
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix");
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(int , t_num, 't',true,"thread num . default[8]");
    END_PARSE_ARGS
    if(! t_num.setted )
    {
        t_num.setted = true ;
        t_num.d.i = 8 ;
    }
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

    {
        BGIQD::MultiThread::MultiThread t_jobs;
        t_jobs.Start(t_num.to_int());
        index  = 0;
        for( auto j : config.keys )
        {
            t_jobs.AddJob([&config, j](){
                    findConnection(config,j,false);
                    }
                    );
            if( config.edge_array[j].id != config.edge_array[j].bal_id )
            {
                unsigned int k = config.edge_array[j].bal_id;
                t_jobs.AddJob([&config, k](){
                        findConnection(config, k,true);
                        }
                        );
            }
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }
    //lger<<BGIQD::LOG::lstart()<<"buildConnection start ... "<<BGIQD::LOG::lend();
    //BGIQD::SOAP2::buildConnection(config);
    lger<<BGIQD::LOG::lstart()<<"linear start ... "<<BGIQD::LOG::lend();
    //BGIQD::SOAP2::LinearConnection(config);
    {
        for( const auto & m : config.keys)
        {
            config.key_array[config.key_map[m]].SetType();
        }
        BGIQD::FREQ::Freq<int> freq;
        BGIQD::FREQ::Freq<int> from;
        BGIQD::FREQ::Freq<int> to;
        for( const auto & m : config.keys )
        {
            if( config.key_array[config.key_map[m]].IsSingle() )
                freq.Touch(0);
            else if ( config.key_array[config.key_map[m]].IsLinear() )
                freq.Touch(2);
            else if ( config.key_array[config.key_map[m]].IsTipTo() )
                freq.Touch(-1);
            else if ( config.key_array[config.key_map[m]].IsTipFrom() )
                freq.Touch(1);
            else 
                freq.Touch(3);
            from.Touch(config.key_array[config.key_map[m]].from.size());
            to.Touch(config.key_array[config.key_map[m]].to.size());
        }
        lger<<BGIQD::LOG::lstart()<<"key type freq"<<'\n'<< freq.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"from freq"<< '\n'<<from.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"to freq"<< '\n'<<to.ToString() << BGIQD::LOG::lend();

        BGIQD::MultiThread::MultiThread t_jobs;
        index = 0;
        t_jobs.Start(t_num.to_int());
        for(size_t i = 1 ; i<= config.keys.size() ; i++)
        {
            t_jobs.AddJob([&config,i](){ linearConnection(config,i); });
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }

    lger<<BGIQD::LOG::lstart()<<"report start ... "<<BGIQD::LOG::lend();
    report(config);
    return 0;
}
