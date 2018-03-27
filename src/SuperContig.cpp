#include "soap2/contigGraph.h"
#include "soap2/loadGraph.h"
#include <iostream>
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"
#include <atomic>

BGIQD::LOG::logger lger;
std::atomic<int> index ;


void findConnection(BGIQD::SOAP2::GlobalConfig & config
        , unsigned int edge_id
        , bool order
        , bool detail = false 
        )
{
    BGIQD::SOAP2::Edge & root = config.edge_array[edge_id];
    unsigned int path_i = edge_id;
    //unsigned int base_i = order ? config.edge_array[edge_id].bal_id : edge_id ;
    std::list<BGIQD::SOAP2::Edge> stack;
    std::map<unsigned int , BGIQD::SOAP2::Edge > history;
    std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > paths;
    std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > mids;
    if ( order )
        config.edge_array[path_i].DepthSearch( config.edge_array , stack,
                history, paths , mids ,config.edge_array[path_i].length , config.connections.at(root.bal_id) );
    else 
        config.edge_array[path_i].DepthSearch( config.edge_array , stack,
                history, paths , mids ,config.edge_array[path_i].length , config.connections.at(root.id) );
    //if( detail )
    //{
    //    for( auto const & i : paths )
    //    {
    //        std::cerr<<"@Path\t"<<edge_id<<'\t'<<i.first<<"\t"<<paths.size()<<std::endl;
    //        for( auto j : i.second )
    //        {
    //            int length = 0 ;
    //            std::stack<BGIQD::SOAP2::Edge> go;
    //            while( !j.empty() )
    //            {
    //                go.push(j.front());
    //                length+=j.top().length;
    //                j.pop();
    //            }
    //            std::cerr<<length;
    //            while(! go.empty() )
    //            {
    //                std::cerr<<go.top().id<<"\t";
    //                go.pop();
    //            }
    //            std::cerr<<"\n";
    //        }
    //    }

    //    for( auto & i : mids)
    //    {
    //        std::cerr<<"@Mids\t"<<edge_id<<'\t'<<i.first<<"\t"<<paths.size()<<std::endl;
    //        for( auto & j : i.second )
    //        {
    //            int length = 0 ;
    //            std::stack<BGIQD::SOAP2::Edge> go;
    //            while( !j.empty() )
    //            {
    //                go.push(j.top());
    //                length+=j.top().length;
    //                j.pop();
    //            }
    //            std::cerr<<length;
    //            while( go.empty() )
    //            {
    //                std::cerr<<go.top().id<<"\t";
    //                go.pop();
    //            }
    //            std::cerr<<"\n";
    //        }
    //    }
    //}
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
        unsigned int to_id_in_path = j.second[0].front().id ;
        //
        //  A1->B1
        //  A2<-B2
        //
        if( !order && to_id == to_id_in_path)
        {
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[path_i]]);
                BGIQD::SOAP2::KeyConn conn{j.first ,0,0};
                conn.SetPostive();
                //conn.length = 
                config.key_array[config.key_map[path_i]].to[j.first]= conn;
            }
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                BGIQD::SOAP2::KeyConn conn{path_i,0,0};
                conn.SetPostive();
                config.key_array[config.key_map[j.first]].from[(path_i)] = conn;
            }
        }
        //
        // A1->B2
        // A2<-B1
        //
        else if( !order && to_id != to_id_in_path )
        {
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[path_i]]);
                BGIQD::SOAP2::KeyConn conn{j.first,0,0};
                config.key_array[config.key_map[path_i]].to[j.first]= conn;
            }
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                BGIQD::SOAP2::KeyConn conn{path_i,0,0};
                config.key_array[config.key_map[j.first]].to[(path_i)] = conn;
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
                BGIQD::SOAP2::KeyConn conn{j.first,0,0};
                config.key_array[config.key_map[root.bal_id]].from[j.first]= conn;
            }
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                BGIQD::SOAP2::KeyConn conn{root.bal_id,0,0};
                config.key_array[config.key_map[j.first]].from[(root.bal_id)] = conn;
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
                BGIQD::SOAP2::KeyConn conn{j.first,0,0};
                conn.SetPostive();
                config.key_array[config.key_map[root.bal_id]].from[j.first]= conn;
            }
            {
                std::lock_guard<std::mutex> lm(config.key_mutex[config.key_map[j.first]]);
                BGIQD::SOAP2::KeyConn conn{root.bal_id,0,0};
                conn.SetPostive();
                config.key_array[config.key_map[j.first]].to[(root.bal_id)] = conn;
            }
        }
    }
    if( index %100 == 0 )
    {
        std::lock_guard<std::mutex> lm(config.contig_mutex);
        lger<<BGIQD::LOG::lstart()<<"process "<<index<<" ..."<<BGIQD::LOG::lend();
    }
}



//
// if 
//  O->A->B && O->B
// delete 
//   O->C
// return the number of deleted conns.
// 
int deleteConns(BGIQD::SOAP2::GlobalConfig &config)
{
    int count = 0;

    auto flush_map = [&config,&count]( std::map<unsigned int , BGIQD::SOAP2::KeyConn> & map , bool order)
    {
        std::vector<BGIQD::SOAP2::KeyConn*> vecs;
        for( auto & i: map )
        {
            vecs.emplace_back(&i.second);
        }
        for( size_t b = 0 ;b< vecs.size() ; b++ )
        {
            for( size_t a = b+1; a< vecs.size() ; a++ )
            {
                if( a == b )
                    continue;
                if( vecs[a]->IsJumpConn() &&  vecs[b]->IsJumpConn() )
                    continue;
                auto & B = config.key_array[config.key_map[vecs[b]->to]];
                auto & A = config.key_array[config.key_map[vecs[a]->to]];
                if( A.IsCircle() || B.IsCircle() )
                    continue;
                bool f11 ,f21 ,f31 ;
                bool f12 ,f22 ,f32 ;
                std::tie(f11,f21,f31) = B.Relationship(vecs[a]->to) ;
                std::tie(f12,f22,f32) = A.Relationship(vecs[b]->to) ;
                if( !f11 || !f12 )
                    continue;
                bool useA1 = false;
                if( vecs[a]->IsPositive() )
                {
                    useA1 = true ;
                }

                if( order )
                {
                    if( useA1 )
                    {
                        //O1->A1 O1->B1 A1->B2
                        //O1->A1 O1->B2 A1->B1
                        //O1->A1 O1->B1 A1<-B2
                        //O1->A1 O1->B2 A1<-B1
                        if ( f32 != vecs[b]->IsPositive() )
                            continue ;
                        if( f22 ) 
                            // O1->A1 O1->B1 A1->B1 --> O1->A1->B1 ;
                            // O1->A1 O1->B2 A1->B2 --> O1->A1->B2 ;
                            // delete O->B
                        {    if ( ! vecs[b]->IsJumpConn() )
                            {
                                vecs[b]->SetJump();
                                count ++;
                            }
                            else {}
                        }
                        else 
                            //O1->A1 O1->B1 A1<-B1 --> O1->B1->A1 ;
                            //O1->A1 O1->B2 A1<-B2 --> O1->B2->A1 ; 
                            //delete O->A
                        {
                            if ( ! vecs[a]->IsJumpConn() )
                            {
                                vecs[a]->SetJump();
                                count ++;
                            }
                            else {}
                        }
                    }
                    else
                    {
                        //O1->A2 O1->B2 A1->B2
                        //O1->A2 O1->B2 A1<-B2
                        //O1->A2 O1->B1 A1->B1
                        //O1->A2 O1->B1 A1<-B1
                        if( f32 == vecs[b]->IsPositive() )
                            continue ;
                        if( f22 ) 
                            // O1->A2 O1->B1 A1->B2 --> O1->B1->A2 ;
                            // O1->A2 O1->B2 A1->B1 --> O1->B2->A2;
                            // delete O->A
                        {
                            if ( ! vecs[a]->IsJumpConn() )
                            {   vecs[a]->SetJump(); count ++; }
                            else {}
                        }
                        else  
                            // O1->A2 O1->B1 A1<-B2 --> O1->A2->B1
                            // O1->A2 O1->B2 A1<-B1 --> O1->A2->B2
                            // delete O->B
                        {
                            if ( ! vecs[b]->IsJumpConn() )
                            {   vecs[b]->SetJump(); count ++; }
                            else {}
                        }
                    }
                }
                else
                {
                    if( useA1)
                    {
                        //O1<-A1 O1<-B1 A1->B2
                        //O1<-A1 O1<-B1 A1<-B2
                        //O1<-A1 O1<-B2 A1->B1
                        //O1<-A1 O1<-B2 A1<-B1
                        if ( f32 != vecs[b]->IsPositive() )
                            continue;
                        if ( f22 )
                        //O1<-A1 O1<-B1 A1->B1 --> A1->B1->O1
                        //O1<-A1 O1<-B2 A1->B2 --> A1->B2->O1
                        //delete O-A
                        {
                            if ( ! vecs[a]->IsJumpConn() )
                            {   vecs[a]->SetJump(); count ++; }
                            else {}
                        }
                        //O1<-A1 O1<-B1 A1<-B1 --> O1<-A1<-B1
                        //O1<-A1 O1<-B2 A1<-B2 --> O1<-A1<-B2
                        //delete O-B
                        else
                        {
                            if ( ! vecs[b]->IsJumpConn() )
                            {   vecs[b]->SetJump(); count ++; }
                            else {}
                        }
                    }
                    else
                    {
                        // O1<-A2 O1<-B1 A1->B1
                        // O1<-A2 O1<-B2 A1->B2
                        // O1<-A2 O1<-B1 A1<-B1
                        // O1<-A2 O1<-B2 A1<-B2
                        if ( f32 == vecs[b]->IsPositive() )
                            continue ;
                        if (! f22 )
                        // O1<-A2 O1<-B1 A1<-B2  --> A2->B1->O1
                        // O1<-A2 O1<-B2 A1<-B1  --> A2->B2->O1
                        // delete O-A
                        {
                            if ( ! vecs[a]->IsJumpConn() )
                            {   vecs[a]->SetJump(); count ++; }
                            else {}
                        }
                        else
                        // O1<-A2 O1<-B1 A1->B2 --> B1->A2->O1
                        // O1<-A2 O2<-B2 A1->B1 --> B2->A2->O1
                        // delete O-B
                        {
                            if ( ! vecs[b]->IsJumpConn() )
                            {   vecs[b]->SetJump(); count ++; }
                            else {}
                        }
                    }
                }
            }
        }
    };

    for( auto i: config.keys )
    {
        auto & curr = config.key_array[config.key_map[i]];
        if( curr.IsCircle() )
            continue;
        if( curr.to.size() > 1 )
        {
            flush_map(curr.to,true);
        }
        if ( curr.from.size() > 1)
        {
            flush_map(curr.from,false);
        }
        index ++ ;
        if( index % 1000 == 0 )
        {
            lger<<BGIQD::LOG::lstart()<<"process "<<index<<" ..."<<BGIQD::LOG::lend();
        }
    }
    return count;
}

int deleteNotBiSupport(BGIQD::SOAP2::GlobalConfig &config)
{
    int count = 0;
    for(auto & i: config.keys)
    {
        auto & curr = config.key_array[config.key_map[i]];
        if( curr.IsCircle() )
            continue;

        for( auto & i : curr.to )
        {
            if( i.second.IsJumpConn() )
                continue;
            const auto & next = config.key_array[config.key_map[i.second.to]];
            bool f1 , f2 , f3 ;
            if( i.second.IsPositive() )
                std::tie(f1,f2,f3) = next.Relationship_nojump( curr.edge_id , false ) ;
            else 
                std::tie(f1,f2,f3) = next.Relationship_nojump( curr.edge_id , true) ;

            if( ! f1 )
            {
                i.second.SetBiNotSuppert();
                count ++;
            }
        }

        for( auto & i : curr.from )
        {
            if( i.second.IsJumpConn() )
                continue;
            const auto & next = config.key_array[config.key_map[i.second.to]];
            bool f1 , f2 , f3 ;
            if( i.second.IsPositive() )
                std::tie(f1,f2,f3) = next.Relationship_nojump( curr.edge_id , true) ;
            else 
                std::tie(f1,f2,f3) = next.Relationship_nojump( curr.edge_id , false) ;
            if( ! f1 )
            {
                i.second.SetBiNotSuppert();
                count ++;
            }
        }
    }
    return count;
}
void linearConnection(BGIQD::SOAP2::GlobalConfig &config , unsigned int key_id)// , BGIQD::MultiThread::MultiThread & queue)
{

    unsigned int i = key_id;
    BGIQD::SOAP2::KeyEdge & curr =  config.key_array[i];
    if( curr.IsMarked() 
            ||curr.IsLinear() || curr.IsCircle()
      )
        return ;
    index ++ ;
    if( curr.IsSingle() )
    {
        curr.Mark();
    }
    else
    {
        BGIQD::SOAP2::ContigRoad path;
        auto extractPath = [&config,&curr,&path]( unsigned int to , bool to_order, bool search_order)
        {
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
            unsigned int next_k = config.key_map[to];
            if( search_order ^ to_order )
            {
                path.real_contig.push_back(config.key_array[next_k].bal_id);
            }
            else
            {
                path.real_contig.push_back(config.key_array[next_k].edge_id);
            }
            bool order = search_order;//
            bool torder = to_order; // 
            if( config.key_array[next_k].IsMarked())
                return ;
            auto  get_nonjump = [] ( const std::map<unsigned int , BGIQD::SOAP2::KeyConn> & map)
            {
                for(const auto & i : map)
                {
                    if( i.second.IsValid() )
                        return std::ref(i.second);
                }
                assert(0);
                return std::ref(map.begin()->second);
            };
            while( (! config.key_array[next_k].IsCircle())&& config.key_array[next_k].IsLinear() && !config.key_array[next_k].IsMarked() )
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
                        const auto & conn  = get_nonjump(config.key_array[next_k].to ) ;
                        next_i = conn.get().to;
                        torder = conn.get().IsPositive();;
                        if( torder )
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].edge_id );
                        }
                        else
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].bal_id);
                        }
                    }
                    //A1->B2
                    else
                    {
                        const auto & conn  = get_nonjump(config.key_array[next_k].from) ;
                        next_i = conn.get().to;
                        torder = conn.get().IsPositive();;
                        if( torder )
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].bal_id);
                        }
                        else
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].edge_id);
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
                        const auto & conn  = get_nonjump(config.key_array[next_k].from) ;
                        next_i = conn.get().to;
                        torder = conn.get().IsPositive();;

                        if( torder )
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].bal_id);
                        }
                        else
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].edge_id);
                        }

                        order = false;
                    }
                    //A1<-B2
                    //A2->B1
                    else
                    {
                        const auto & conn  = get_nonjump(config.key_array[next_k].to ) ;
                        next_i = conn.get().to;
                        torder = conn.get().IsPositive();;
                        order = true;

                        if( torder )
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].edge_id);
                        }
                        else
                        {
                            path.real_contig.push_back(config.key_array[config.key_map[next_i]].bal_id);
                        }
                    }
                }
                next_k = config.key_map[next_i];
                path.contig.push_back(next_i) ;
            }
            if( config.key_array[next_k].IsMarked() )
                return ;
            // line->next_k->
            // <-next_k<-line
            path.tailin = false;
            if(! config.key_array[next_k].IsCircle() )
            {

                if( ( torder && order ) || ( !torder && !order ))
                {
                    if( config.key_array[next_k].from_size == 1 )
                    {
                        path.tailin = true;
                    }
                }
                else
                {
                    if( config.key_array[next_k].to_size == 1 )
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
                std::lock_guard<std::mutex> lm(config.contig_mutex);
                config.contigs.push_back(path);
            }
        };
        for( auto next : curr.to )
        {
            if(! next.second.IsValid() )
                continue;
            extractPath(next.first,next.second.IsPositive(),true);
        }
        for(auto next : curr.from)
        {
            if( ! next.second.IsValid() )
                continue;
            extractPath(next.first,next.second.IsPositive(),false);
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


    //std::cout<<"--- Paths start  ----"<<std::endl;
    for(const auto & i : config.contigs)
    {
        freq.Touch(i.length);
        if( i.downstream )
            std::cout<<'>';
        else
            std::cout<<'<';

        if( i.headin )
            std::cout<<'[';
        else
            std::cout<<'(';
        std::cout<<*i.contig.begin()<<'\t'<<*i.contig.rbegin();
        if( i.tailin )
            std::cout<<']';
        else
            std::cout<<')';
        std::cout<<'\t'<<i.length<<'\t';

        for( auto j : i.contig)
            std::cout<<j<<'\t';
        std::cout<<std::endl;
    }

    auto fout= BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(config.contigroad);
    for(const auto & i : config.contigs)
    {
        if( i.length < 2 )
            continue;
        (*fout)<<i.length<<'\t';
        int start = 1 ;
        if( i.headin )
        {
            start --;
        }
        for( int j =0 ; j<i.length ; j++ )
        {
            (*fout)<<i.real_contig[start+j]<<"\t";
        }
        (*fout)<<std::endl;
    }
    delete fout;

    lger<<BGIQD::LOG::lstart()<<"clusterNum "<<config.clusterNum<<BGIQD::LOG::lend();
    lger<<BGIQD::LOG::lstart()<<"pathNum "<<config.contigs.size()<<BGIQD::LOG::lend();
    lger<<BGIQD::LOG::lstart()<<"linear length freq\n"<<freq.ToString()<<BGIQD::LOG::lend();
}

int main(int argc , char **argv)
{
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , lger);
    BGIQD::LOG::timer t(lger,"SuperContig");
    START_PARSE_ARGS
        DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix");
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(int , t_num, 't',true,"thread num . default[8]");
    DEFINE_ARG_DETAIL(bool , super, 's',true,"super contig ? default false");
    DEFINE_ARG_DETAIL(bool , detail, 'd',true,"print detail ? default false");
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
    config.contigroad =  prefix.to_string() +".contigroad";

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
            t_jobs.AddJob([&config, j,&detail](){
                    findConnection(config,j,false , detail.to_bool());
                    }
                    );
            if( config.edge_array[j].id != config.edge_array[j].bal_id )
            {
                unsigned int k = config.edge_array[j].bal_id;
                t_jobs.AddJob([&config, k,&detail](){
                        findConnection(config, k,true ,detail.to_bool());
                        }
                        );
            }
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }
    //lger<<BGIQD::LOG::lstart()<<"buildConnection start ... "<<BGIQD::LOG::lend();
    //BGIQD::SOAP2::buildConnection(config);
    for( const auto & m : config.keys)
    {
        config.key_array[config.key_map[m]].CheckCircle();
    }
    lger<<BGIQD::LOG::lstart()<<"deleteConn start ... "<<BGIQD::LOG::lend();
    {
        index = 0 ;
        //while(1)
        //{
        int dd1 = deleteNotBiSupport(config);
        lger<<BGIQD::LOG::lstart()<<"deleteNotBiSupport "<<dd1<<" conn"<<BGIQD::LOG::lend();
        int d = deleteConns(config);
        lger<<BGIQD::LOG::lstart()<<"deleteConn delete  "<<d<<" conn"<<BGIQD::LOG::lend();

        //    if( d == 0 )
        //        break;
        //}
        int dd = deleteNotBiSupport(config);
        lger<<BGIQD::LOG::lstart()<<"deleteNotBiSupport "<<dd<<" conn"<<BGIQD::LOG::lend();

    }
    lger<<BGIQD::LOG::lstart()<<"freq report... "<<BGIQD::LOG::lend();
    {
        for( const auto & m : config.keys)
        {
            config.key_array[config.key_map[m]].SetType();
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

        for( const auto & m : config.keys )
        {

            auto & curr = config.key_array[config.key_map[m]];
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


            if( config.key_array[config.key_map[m]].IsSingle() )
                freq.Touch("Single");
            else if ( config.key_array[config.key_map[m]].IsLinear() )
                freq.Touch("Linear");
            else if ( config.key_array[config.key_map[m]].IsTipTo() && config.key_array[config.key_map[m]].from.size() ==1 )
                freq.Touch("Tipto 1");
            else if ( config.key_array[config.key_map[m]].IsTipTo() && config.key_array[config.key_map[m]].from.size() > 1 )
                freq.Touch("Tipto >1");
            else if ( config.key_array[config.key_map[m]].IsTipFrom() && config.key_array[config.key_map[m]].to.size() ==1 )
                freq.Touch("Tipfrom 1");
            else if ( config.key_array[config.key_map[m]].IsTipFrom() && config.key_array[config.key_map[m]].to.size() > 1 )
                freq.Touch("Tipfrom >1");
            else
                freq.Touch("Multi");
            from.Touch(config.key_array[config.key_map[m]].from_size);
            to.Touch(config.key_array[config.key_map[m]].to_size);
            total.Touch(config.key_array[config.key_map[m]].total_size);
            del.Touch(config.key_array[config.key_map[m]].jump_conn);
            base_from.Touch( config.key_array[config.key_map[m]].from.size());
            base_to.Touch( config.key_array[config.key_map[m]].to.size());
            in_circle.Touch(config.key_array[config.key_map[m]].IsCircle());
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
    if ( ! super.to_bool() )
        return 0;
    lger<<BGIQD::LOG::lstart()<<"linear start ... "<<BGIQD::LOG::lend();
    //BGIQD::SOAP2::LinearConnection(config);
    {
        BGIQD::MultiThread::MultiThread t_jobs;
        index = 0;
        //TODO : make it thread safe
        t_jobs.Start(1);
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
