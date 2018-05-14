#include "stLFR/barcodeOnContig.h"
#include "stLFR/ContigCluster.h"
#include "stLFR/LineGroup.h"
#include "soap2/contigGraph.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_reader.h"
#include "common/freq/freq.h"
#include "soap2/contigGraphSPF.h"

#include <set>
#include <queue>

typedef BGIQD::GRAPH::EdgeIterator<BGIQD::SOAP2::GraphEA_Access> EdgeItr;

typedef BGIQD::GRAPH::SPFSearch<
BGIQD::SOAP2::GraphEA_Access,
    EdgeItr,
    BGIQD::SOAP2::SFPEnder
    > Searcher;

typedef Searcher::SPFNode SNode;

struct GlobalConfig
{
    BGIQD::stLFR::GraphEA_withBarcode graph_eab;
    BGIQD::stLFR::ContigCluster clusters;
    BGIQD::stLFR::ContigRoads roads;
    BGIQD::LOG::logger lger;


    int K;
    std::string updateEdge;
    std::string arc;
    std::string cluster;
    std::string road;
    std::string read2contig;

    void LoadKeyInfo()
    {
        for(const auto & i : clusters.connections)
        {
            graph_eab.graph_ea.edge_array[i.first].SetKey();
        }
    }
    BGIQD::FREQ::Freq<int> path_num_freq;
    BGIQD::FREQ::Freq<int> circle_run;
    BGIQD::FREQ::Freq<std::string> road_fill_freq;

} config;

struct SearchResult
{
    //std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > paths;
    //std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > mids;
    bool rearch;

    bool head_tail ; // depth search from head or tail
    bool downstream;

    enum Status
    {
        NotRearch = 0,
        A1B1_B2A2 = 1 ,
        A1B2_B1A2 = 2 ,
        B1A1_A2B2 = 3 ,
        B2A1_A2B1 = 4 ,
    } status;

    unsigned int true_from ;
    unsigned int true_to ;

    Searcher searcher;
};

SearchResult SearchAllPath(unsigned int from  , unsigned int to, int max){

    SearchResult ret;
    ret.head_tail = true ;
    ret.downstream = true ;
    ret.status = SearchResult::NotRearch;
    BGIQD::SOAP2::Edge * edge_array = config.graph_eab.graph_ea.edge_array;

    unsigned int head_id = from;
    unsigned int tail_id = to;
    unsigned int search_id = head_id;
    unsigned int target_id = tail_id;

    std::map<unsigned int ,float> neibs ;
    auto & node = edge_array[target_id];
    unsigned int to_key_id ;
    if( node.IsKey() )
        to_key_id = target_id ;
    else
        to_key_id = node.bal_id ;

    neibs[to_key_id] = 0.5f;

    auto key = [neibs, edge_array] (unsigned int id)
    {
        auto & node = edge_array[id] ;
        auto & node_bal = edge_array[node.bal_id];
        if( ! node.IsKey() && !node_bal.IsKey() )
        {
            return BGIQD::SOAP2::NodeType::Normal ;
        }
        if ( node.IsKey() )
        {
            if( neibs.find(id) == neibs.end() )
                return BGIQD::SOAP2::NodeType::Key_Unknow ;
            else
                return BGIQD::SOAP2::NodeType::Key_Neibs ;
        }
        else
        {
            if( neibs.find(node.bal_id) == neibs.end() )
                return BGIQD::SOAP2::NodeType::RC_Key_Unknow ;
            else
                return BGIQD::SOAP2::NodeType::RC_Key_Neibs ;
        }
    };

    ret.searcher.accesser.base = &config.graph_eab.graph_ea;
    ret.searcher.ender.Init( key , max);

    ret.searcher.DoSPFSearch(search_id);

    ret.true_from = from ;
    ret.true_to = to ;

    auto & ender = ret.searcher.ender ;
    auto itr = ender.founder.find(to_key_id);
    if( itr == ender.founder.end() )
    {
        assert(0);
        ret.status = SearchResult::NotRearch;
        return ret;
    }
    auto & tos = itr->second;
    if( node.IsKey() &&  ! tos.base )
    {
        assert(0);
        ret.status = SearchResult::NotRearch;
        return ret;
    }
    else if ( ! node.IsKey() && ! tos.bal )
    {
        assert(0);
        ret.status = SearchResult::NotRearch;
        return ret;
    }
    else if ( node.IsKey() && tos.base )
    {
        ret.status = SearchResult::A1B1_B2A2 ;

    }
    else if ( ! node.IsKey() && tos.bal ) 
    {
        ret.status = SearchResult::A1B2_B1A2 ;

    }
    else
    {
        assert(0);
        ret.status = SearchResult::NotRearch;
        return ret;
    }
    return ret ;
}

void FindCorrectPath(unsigned int from , unsigned int to, 
        const SearchResult & result , BGIQD::stLFR::P2PGraph & p2pgrapg
        )
{
    p2pgrapg.Init(from,to);

    std::set<unsigned int> history;
    std::queue<SNode> nexts;

    try{
        auto & end = result.searcher.fib_nodes.at(to);
        history.insert(to) ;
        nexts.push(end);
    }
    catch( ... ) 
    {
        assert(0);
    }

    try{
        while(! nexts.empty() )
        {
            auto & a_to = nexts.front() ;
            if( a_to.Base.value != from )
            {
                p2pgrapg.AddFromTo( a_to.prev , a_to.Base.value ) ;//, a_to.prev );
                if( history.find( a_to.prev ) == history.end() )
                {
                    history.insert( a_to.prev ) ;
                    nexts.push( result.searcher.fib_nodes.at(a_to.prev));
                }
            }
            for( auto from : a_to.other_from )
            {
                p2pgrapg.AddFromTo(from , a_to.Base.value ) ;//,from );
                if( history.find( from ) == history.end() )
                {
                    history.insert( from ) ;
                    nexts.push( result.searcher.fib_nodes.at(from));
                }
            }
            nexts.pop();
        }
    }
    catch( ... )
    {
        assert(0);
    }
    p2pgrapg.GeneratePath();
}

bool AppendPath( const  BGIQD::stLFR::P2PGraph & p2pgrapg , const SearchResult & result, BGIQD::stLFR::ContigRoad & road)
{
    if( ! road.needMerge() )
        return false;
    if( result.status== SearchResult::NotRearch )
        return false;
    if( p2pgrapg.path_num < 1 )
        return false;

    if( result.head_tail )
    {
        road.contig_path.insert(road.contig_path.end() 
                , p2pgrapg.final_path.begin()
                , p2pgrapg.final_path.end() );
    }
    else
    {
        road.contig_path.insert(road.contig_path.end() 
                , p2pgrapg.final_path.rbegin()
                , p2pgrapg.final_path.rend() );
    }

    road.circle_runs.push_back(p2pgrapg.final_circled);

    return true;
}

void FillContigRoad( BGIQD::stLFR::ContigRoad & road, int max , float ecov , bool circle_solve)
{
    static std::mutex write_mutex;
    static std::mutex path_num_mutex;
    if ( ! road.needMerge() )
        return ;

    road.status = BGIQD::stLFR::ContigRoad::FillStatus::None ;
    road.fill_num = 0;
    bool line_down = false;
    for(int i=1  ; i < road.linear_length ; i++ )
    {
        auto start = road.getLinearStep(i);
        auto ret = SearchAllPath(start.first , start.second,max);
        if( ret.status == SearchResult::NotRearch )
        {
            if( road.status  == BGIQD::stLFR::ContigRoad::FillStatus::None )
            {
                continue ;
            }
            else
            {
                break;
            }
        }
        else
        {
            bool curr_line_down = false ;

            if ( ret.status == SearchResult::A1B1_B2A2 ||  ret.status == SearchResult::A1B2_B1A2  )
            {
                curr_line_down = true ;
            }

            if( road.status  == BGIQD::stLFR::ContigRoad::FillStatus::None )
            {
                line_down = curr_line_down ;
            }
            else if ( line_down != curr_line_down )
            {
                road.status = BGIQD::stLFR::ContigRoad::FillStatus::Conflict ;
                break ;
            }

            /*BGIQD::stLFR::P2PGraph p2pgrapg;

            p2pgrapg.base_graph = &config.graph_eab;
            p2pgrapg.deal_circle  = circle_solve ;
            p2pgrapg.ecov = ecov ;
            p2pgrapg.K =config.K ;
            FindCorrectPath(ret.true_from , ret.true_to , ret , p2pgrapg );
            // check if allpath find a correct path ?
            {
                std::lock_guard<std::mutex> l(path_num_mutex);
                config.path_num_freq.Touch( p2pgrapg.path_num);
            }

            if( p2pgrapg.path_num < 1 )
            {
                if( road.status  == BGIQD::stLFR::ContigRoad::FillStatus::None )
                {
                    continue ;
                }
                else
                {
                    break;
                }
            }
            */
            auto extractShortestPath = [] (const Searcher & searcher , unsigned int from , unsigned int to)
            {
                unsigned int curr = to ;
                std::stack<unsigned int> path;
                while(curr != from )
                {
                    path.push(curr);
                    try {
                        const auto & node = searcher.fib_nodes.at(curr);
                        curr = node.prev ;
                    }
                    catch ( ... ) 
                    {
                        assert(0);
                    }
                };
                std::vector<unsigned int> ret ;
                while(!  path.empty() )
                {
                    ret.push_back(path.top());
                    path.pop() ;
                }
                return ret ;
            };

            if( road.fill_num== 0 )
            {
                if( ret.head_tail )
                    road.contig_path.push_back( ret.true_from );
                else 
                    road.contig_path.push_back( ret.true_to );
            }

            auto retp = extractShortestPath( ret.searcher , ret.true_from , ret.true_to );
            road.contig_path.insert(road.contig_path.end() 
                    , retp.begin()
                    , retp.end() );
            //if( ! AppendPath( p2pgrapg ,ret,road ) )
            //{
            //    road.status = BGIQD::stLFR::ContigRoad::FillStatus::Conflict ;
            //    break;
            //}
            //if( ret.head_tail )
            //    road.contig_path.push_back( ret.true_to );
            //else 
            //    road.contig_path.push_back( ret.true_from );
            road.status = BGIQD::stLFR::ContigRoad::FillStatus::PartSucc ;
            road.fill_num ++ ;
        }
    }

    if ( road.fill_num == road.linear_length - 1 )
    {
        road.status = BGIQD::stLFR::ContigRoad::FillStatus::Complete;
    }
    /*{
        std::lock_guard<std::mutex> l(write_mutex);
        if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Conflict)
        {
            config.road_fill_freq.Touch("Conflict");
            return;
        }
        if( road.status == BGIQD::stLFR::ContigRoad::FillStatus::None )
        {
            config.road_fill_freq.Touch("None");
            return;
        }
        if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Complete)
        {
            config.road_fill_freq.Touch("Complete");
        }
        if ( road.status == BGIQD::stLFR::ContigRoad::FillStatus::PartSucc)
        {
            config.road_fill_freq.Touch("PartSucc");
        }
        for( const auto i : road.contig_path )
        {
            std::cout<<i<<'\t';
        }
        std::cout<<std::endl;
    }*/
}

void report()
{
    int filled = 0 ;
    for( const auto &road : config.roads.roads)
    {
        if( ! road.needMerge() )
            continue;
        if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Conflict)
        {
            config.road_fill_freq.Touch("Conflict");
            continue;
        }
        if( road.status == BGIQD::stLFR::ContigRoad::FillStatus::None )
        {
            config.road_fill_freq.Touch("None");
            continue;
        }
        if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Complete)
        {
            config.road_fill_freq.Touch("Complete");
        }
        if ( road.status == BGIQD::stLFR::ContigRoad::FillStatus::PartSucc)
        {
            config.road_fill_freq.Touch("PartSucc");
        }
        bool has_circle = false ;
        for( auto i : road.circle_runs )
        {
            config.circle_run.Touch(i);
            if( i > 0 )
                has_circle = true ;
        }
        if( has_circle )
        {
            std::cout<<"c\t";
        }
        else
        {
            std::cout<<"l\t";
        }
        for( const auto i : road.contig_path )
        {
            std::cout<<i<<'\t';
        }
        filled += (road.fill_num + 1 );
        std::cout<<std::endl;
    }
    config.lger<<BGIQD::LOG::lstart()<<"road fill use seed contig \n"<<filled<<BGIQD::LOG::lend();
}

int  main(int argc, char **argv)
{
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , config.lger);
    BGIQD::LOG::timer t(config.lger,"SuperContig");
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , prefix, 'o',false,"prefix \n \
                                            need    xxx.Arc\n\
                                                    xxx.updated.edge\n\
                                                    xxx.cluster\n\
                                                    xxx.contigroad\n\
                                                    xxx.read2contig");
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(int , t_num, 't',true,"thread num . default[8]");
    DEFINE_ARG_DETAIL(float , Ecov, 'e',false,"Ecov of contigs");
    DEFINE_ARG_DETAIL(bool, circle, 'c',true,"solve circle by cov ? default not");
    DEFINE_ARG_DETAIL(int, searchDepth, 'l',true,"search depth (bp) default 10000");
    END_PARSE_ARGS

    config.K = kvalue.to_int();
    config.arc = prefix.to_string() +".Arc";
    config.updateEdge = prefix.to_string() +".updated.edge";
    config.cluster= prefix.to_string() +".cluster";
    config.road = prefix.to_string() +".contigroad";
    config.read2contig = prefix.to_string() +".read2contig";

    if( ! t_num.setted )
    {
        t_num.setted = true ;
        t_num.d.i = 8 ;
    }
    config.lger<<BGIQD::LOG::lstart()<<"parse args end ... "<<BGIQD::LOG::lend();
    //step1 .loading ...
    config.graph_eab.graph_ea.LoadEdge(config.updateEdge,config.K);
    config.lger<<BGIQD::LOG::lstart()<<"load updateEdg end ... "<<BGIQD::LOG::lend();
    config.graph_eab.graph_ea.LoadArc(config.arc);
    config.lger<<BGIQD::LOG::lstart()<<"load arc  end ... "<<BGIQD::LOG::lend();
    config.graph_eab.LoadBarcodeOnConfig(config.read2contig);
    config.lger<<BGIQD::LOG::lstart()<<"load read2contig end ... "<<BGIQD::LOG::lend();
    config.clusters.loadCluster(config.cluster);
    config.lger<<BGIQD::LOG::lstart()<<"load cluster end ... "<<BGIQD::LOG::lend();
    config.roads.LoadRoads(config.road);
    config.lger<<BGIQD::LOG::lstart()<<"load road end ... "<<BGIQD::LOG::lend();
    config.LoadKeyInfo();
    config.lger<<BGIQD::LOG::lstart()<<"load key end ... "<<BGIQD::LOG::lend();

    //step2 fill road ... 
    {
        BGIQD::MultiThread::MultiThread t_jobs;
        t_jobs.Start(t_num.to_int());

        int max = searchDepth.to_int();
        float ecov = Ecov.to_float() ;
        bool circle_solve = circle.to_bool() ;
        for(int i= 0 ; i<(int)config.roads.roads.size(); i++)
        {
            t_jobs.AddJob([i, max ,ecov ,circle_solve](){ FillContigRoad(std::ref(config.roads.roads[i]),max,ecov,circle_solve); });
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }

    config.lger<<BGIQD::LOG::lstart()<<"fill contig road end ... "<<BGIQD::LOG::lend();

    report();

    config.lger<<BGIQD::LOG::lstart()<<"report end ... "<<BGIQD::LOG::lend();
    config.lger<<BGIQD::LOG::lstart()<<"all path freq \n"<<config.path_num_freq.ToString()<<BGIQD::LOG::lend();
    config.lger<<BGIQD::LOG::lstart()<<"road fill freq \n"<<config.road_fill_freq.ToString()<<BGIQD::LOG::lend();
    config.lger<<BGIQD::LOG::lstart()<<"circled  freq \n"<<config.circle_run.ToString()<<BGIQD::LOG::lend();

    /*for( const auto & road :  config.roads.roads )
    {
        if(road.needMerge()
                && ( road.status != BGIQD::stLFR::ContigRoad::FillStatus::Conflict 
                    || road.status != BGIQD::stLFR::ContigRoad::FillStatus::None )) 
        {
            for( const auto i : road.contig_path )
            {
                std::cout<<i<<'\t';
            }
            std::cout<<std::endl;
        }
    }*/
    //step2 print road ... 
}
