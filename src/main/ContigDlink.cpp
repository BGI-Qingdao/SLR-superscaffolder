#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"
#include "common/error/Error.h"

#include <atomic>
#include <iostream>
#include <algorithm>
#include <set>
#include <mutex>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>

#include "soap2/contigGraphDepth.h"
#include "soap2/contigGraphSPF.h"
#include "soap2/contigGraph.h"
#include "soap2/fileName.h"
#include "soap2/contigGraph.h"

#include "common/files/file_reader.h"
#include "common/error/Error.h"

struct AppConfig
{
    BGIQD::SOAP2::GraphEA  graph_ea;
    BGIQD::SOAP2::Edge * edge_array;

    // cluster(key) data
    unsigned int clusterNum;
    std::set<unsigned int> keys;
    std::map<unsigned int , std::map<unsigned int,float > > connections;
    long long connectionNum;

    // edge id --> key id
    std::map<unsigned int , unsigned int> key_map;
    BGIQD::SOAP2::KeyEdge * key_array;
    std::mutex * key_mutex;

    // super conitgs
    std::mutex contig_mutex;
    std::vector<BGIQD::SOAP2::ContigRoad> contigs ;//std::vector<unsigned int> > contigs;

    void loadUpdateEdge( )
    {
        graph_ea.LoadEdge(fName.updatedEdge(),K);
        edge_array = graph_ea.edge_array;
    }
    void loadArc()
    {
        graph_ea.LoadArc(fName.Arc());
    }

    void loadCluster()
    {
        std::string line;
        unsigned int contigId;
        unsigned int to;
        float cov;
        connectionNum= 0;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.cluster());
        if( in == NULL )
            FATAL( " open prefix.cluser for read failed !!! " );
        // load connection 
        while(!std::getline(*in,line).eof())
        {
            std::istringstream ist(line);
            ist>>contigId;
            //config.connections[contigId][contigId] = 1.0f;

            keys.insert(contigId);
            while(! ist.eof() )
            {
                ist>>to>>cov;
                connections[contigId][to] = cov;
                connections[to][contigId] = cov ;
                //config.keys.insert(to);
            }
        }
        delete in ;

        // init keys
        clusterNum = keys.size();
        key_array =static_cast<BGIQD::SOAP2::KeyEdge*>( calloc( sizeof(BGIQD::SOAP2::KeyEdge) , clusterNum +1));
        key_mutex = new std::mutex[clusterNum+1];
        unsigned int index = 1;
        for( const auto & i : keys)
        {
            edge_array[i].SetKey();
            key_array[index] = BGIQD::SOAP2::KeyEdge();
            key_array[index].edge_id = i;
            key_array[index].bal_id = edge_array[i].bal_id;
            key_array[index].id = index;
            key_map[i] = index ;
            index++;
        }
    }
    //void buildConnection(GlobalConfig&);
    //void LinearConnection(GlobalConfig &);
    BGIQD::LOG::logger lger;
    std::atomic<int> index ;
    BGIQD::SOAP2::FileNames fName;
    int K;
    void Init(const std::string & prefix , int k )
    {
        fName.Init(prefix);
        K = k ;
        BGIQD::LOG::logfilter::singleton().get("ContigDlink",BGIQD::LOG::loglevel::INFO , lger);
        index = 0;
        lger<<BGIQD::LOG::lstart()<<"loadUpdateEdge start ... "<<BGIQD::LOG::lend();
        loadUpdateEdge();
        lger<<BGIQD::LOG::lstart()<<"loadArc start ... "<<BGIQD::LOG::lend();
        loadArc();
        lger<<BGIQD::LOG::lstart()<<"loadCluster start ... "<<BGIQD::LOG::lend();
        loadCluster();
    }

    void findConnection(unsigned int edge_id
            , bool is_bal 
            , int max_length
            , int K
            )
    {
        BGIQD::SOAP2::Edge & root = edge_array[edge_id];

        unsigned int real_from = edge_id ;

        unsigned int key_from = is_bal ? root.bal_id : edge_id ;

        std::map<unsigned int , float> * neibs;

        neibs = & connections.at(key_from);

        auto key = [neibs,this] (unsigned int id)
        {
            auto & node = edge_array[id] ;
            auto & node_bal = edge_array[node.bal_id];
            if( ! node.IsKey() && !node_bal.IsKey() )
            {
                return BGIQD::SOAP2::NodeType::Normal ;
            }
            if ( node.IsKey() )
            {
                if( neibs->find(id) == neibs->end() )
                    return BGIQD::SOAP2::NodeType::Key_Unknow ;
                else
                    return BGIQD::SOAP2::NodeType::Key_Neibs ;
            }
            else
            {
                if( neibs->find(node.bal_id) == neibs->end() )
                    return BGIQD::SOAP2::NodeType::RC_Key_Unknow ;
                else
                    return BGIQD::SOAP2::NodeType::RC_Key_Neibs ;
            }
        };

        typedef BGIQD::GRAPH::EdgeIterator<BGIQD::SOAP2::GraphEA_Access> EdgeItr;

        typedef BGIQD::GRAPH::SPFSearch<
            BGIQD::SOAP2::GraphEA_Access,
            EdgeItr,
            BGIQD::SOAP2::SFPEnder
                > Searcher;
        Searcher searcher;
        searcher.accesser.base = &graph_ea ;
        searcher.accesser.K = K ;
        searcher.ender.Init( key , max_length );

        searcher.DoSPFSearch(real_from);
        index ++ ;

        for ( auto & j : searcher.ender.founder )
        {
            unsigned int key_to = j.first ;

            float sim = neibs->at(key_to);

            auto & curr_from = key_array[key_map.at(key_from) ];

            auto & curr_to = key_array[key_map[key_to]];

            //
            //  A1->B1 * +
            //  A2<-B2
            //

            if ( ! is_bal && j.second.base )
            {

                unsigned int real_to = key_to ;

                int length  = searcher.fib_nodes.at(real_to).Base.key + K 
                    - searcher.accesser.AccessNode(real_to).length;

                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_from]]);
                    BGIQD::SOAP2::KeyConn conn{key_to,length,0,sim};
                    conn.SetPostive();
                    curr_from.to[key_to]= conn;
                }
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_to]]);
                    BGIQD::SOAP2::KeyConn conn{key_from,length,0,sim};
                    conn.SetPostive();
                    curr_to.from[key_from] = conn;
                }
            }

            //
            // A1->B2 * +
            // A2<-B1 +
            //

            if( ! is_bal && j.second.bal)
            {
                unsigned int real_to =curr_to.bal_id;
                unsigned int B1_to = curr_from.bal_id ;
                int length  = searcher.fib_nodes.at(real_to).Base.key + K // nodes.at(curr_to.bal_id).path_length
                    - searcher.accesser.AccessNode(real_to).length;
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_from]]);
                    BGIQD::SOAP2::KeyConn conn{key_to,length,0, sim};
                    curr_from.to[real_to]= conn;
                }
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_to]]);
                    BGIQD::SOAP2::KeyConn conn{key_from,length,0, sim};
                    curr_to.to[B1_to] = conn;
                }
            }
            //
            // A2->B1 * +
            // A1<-B2 +
            //
            if ( is_bal && j.second.base )
            {
                unsigned int real_to = key_to ;
                unsigned int A1_to = curr_to.bal_id;
                int length  = searcher.fib_nodes.at(real_to).Base.key + K
                    - searcher.accesser.AccessNode(real_to).length;
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_from]]);
                    BGIQD::SOAP2::KeyConn conn{key_to,length,0, sim};
                    curr_from.from[A1_to]= conn;
                }
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_to]]);
                    BGIQD::SOAP2::KeyConn conn{key_from,length,0, sim};
                    curr_to.from[real_from] = conn;
                }
            }
            //
            // A2->B2 *
            // A1<-B1 +
            //
            if ( is_bal && j.second.bal )
            {
                int real_to =curr_to.bal_id;

                int length  = searcher.fib_nodes.at(real_to).Base.key + K 
                    - searcher.accesser.AccessNode(real_to).length;
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_from]]);
                    BGIQD::SOAP2::KeyConn conn{key_to,length,0, sim};
                    conn.SetPostive();
                    curr_from.from[key_to]= conn;
                }
                {
                    std::lock_guard<std::mutex> lm(key_mutex[key_map[key_to]]);
                    BGIQD::SOAP2::KeyConn conn{key_from,length,0 , sim};
                    conn.SetPostive();
                    curr_to.to[key_from] = conn;
                }
            }
        }
    }

    void LogFreq()
    {

        for( const auto & m : keys)
        {
            key_array[key_map[m]].CheckCircle();
        }

        for( const auto & m : keys)
        {
            key_array[key_map[m]].SetType();
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

        for( const auto & m : keys )
        {

            auto & curr = key_array[key_map[m]];
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


            if( key_array[key_map[m]].IsSingle() )
                freq.Touch("Single");
            else if ( key_array[key_map[m]].IsLinear() )
                freq.Touch("Linear");
            else if ( key_array[key_map[m]].IsTipTo() && key_array[key_map[m]].to_size ==1 )
                freq.Touch("Tipto 1");
            else if ( key_array[key_map[m]].IsTipTo() && key_array[key_map[m]].to_size > 1 )
                freq.Touch("Tipto >1");
            else if ( key_array[key_map[m]].IsTipFrom() && key_array[key_map[m]].from_size ==1 )
                freq.Touch("Tipfrom 1");
            else if ( key_array[key_map[m]].IsTipFrom() && key_array[key_map[m]].from_size > 1 )
                freq.Touch("Tipfrom >1");
            else
                freq.Touch("Multi");
            from.Touch(key_array[key_map[m]].from_size);
            to.Touch(key_array[key_map[m]].to_size);
            total.Touch(key_array[key_map[m]].total_size);
            del.Touch(key_array[key_map[m]].jump_conn);
            base_from.Touch( key_array[key_map[m]].from.size());
            base_to.Touch( key_array[key_map[m]].to.size());
            in_circle.Touch(key_array[key_map[m]].IsCircle());
        }

        lger<<BGIQD::LOG::lstart()<<"base key type freq"<<'\n'<<basefreq.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"key type freq"<<'\n'<< freq.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"from freq"<< '\n'<<from.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"to freq"<< '\n'<<to.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"total freq"<< '\n'<<total.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"delete freq"<< '\n'<<del.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"basefrom freq"<< '\n'<<base_from.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"baseto freq"<< '\n'<<base_to.ToString() << BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"incircle freq"<< '\n'<<in_circle.ToString() << BGIQD::LOG::lend();
    }

    void PrintContigDlinkGraph()
    {
        auto deg= BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.connInfo());
        if( deg == NULL )
            FATAL( "open prefix.connInfo for write failed !!! ");
        for( const auto & m : keys )
        {

            auto & curr = key_array[key_map[m]];
            if( curr.from.size() > 0 )
            {
                (*deg)<<curr.edge_id<<"\t-\t";
                for( auto i: curr.from )
                {
                    (*deg)<<i.second.to
                        <<":"<<i.second.length
                        <<":"<<i.second.sim
                        <<":"<<(i.second.IsPositive() ? '+' :'-')<<"\t";
                }
                *(deg)<<std::endl;
            }
            if( curr.to.size() > 0 )
            {
                (*deg)<<curr.edge_id<<"\t+\t";
                for( auto i: curr.to)
                {
                    (*deg)<<i.second.to
                        <<":"<<i.second.length
                        <<":"<<i.second.sim
                        <<":"<<(i.second.IsPositive() ? '+' :'-')<<"\t";
                }
                *(deg)<<std::endl;
            }
        }
        delete deg;
    }

} config;



int main(int argc , char **argv)
{
    START_PARSE_ARGS;
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix. Input xxx.cluster && xxx.Arc && xxx.update.edge ; Ouput xxx.connInfo" );
    DEFINE_ARG_REQUIRED(int , kvalue,"K value");
    DEFINE_ARG_OPTIONAL(int , thread, "thread num .","8");
    DEFINE_ARG_OPTIONAL(int, searchDepth,"search depth (bp) ","10000");
    END_PARSE_ARGS;


    config.Init(prefix.to_string() , kvalue.to_int());

    config.lger<<BGIQD::LOG::lstart()<<"parse args end ... "<<BGIQD::LOG::lend();

    BGIQD::LOG::timer t(config.lger,"ContigDlink");

    config.lger<<BGIQD::LOG::lstart()<<"buildConnection start ... "<<BGIQD::LOG::lend();
    {
        BGIQD::MultiThread::MultiThread t_jobs;
        t_jobs.Start(thread.to_int());
        int searchDepth_value = searchDepth.to_int();
        for( auto j : config.keys )
        {
            t_jobs.AddJob([j, searchDepth_value](){
                    config.findConnection(j ,false  , searchDepth_value, config.K);
                    }
                    );
            if( config.edge_array[j].id != config.edge_array[j].bal_id )
            {
                unsigned int k = config.edge_array[j].bal_id;
                t_jobs.AddJob([ k, searchDepth_value](){
                        config.findConnection(k ,true , searchDepth_value, config.K);
                        }
                        );
            }
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }

    config.PrintContigDlinkGraph();

    config.LogFreq();

    return 0;
}
