#include "stLFR/barcodeOnContig.h"
#include "stLFR/ContigCluster.h"
#include "stLFR/LineGroup.h"

#include "soap2/contigGraph.h"
#include "soap2/contigGraphSPF.h"
#include "soap2/fileName.h"

#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/freq/freq.h"
#include "common/error/Error.h"
#include <set>
#include <queue>

typedef BGIQD::GRAPH::EdgeIterator<BGIQD::SOAP2::GraphEA_Access> EdgeItr;

typedef BGIQD::GRAPH::SPFSearch<
BGIQD::SOAP2::GraphEA_Access,
    EdgeItr,
    BGIQD::SOAP2::SFPEnder
    > Searcher;

typedef Searcher::SPFNode SNode;

struct SearchResult
{
    bool succ ;

    unsigned int true_from ;
    unsigned int true_to ;

    Searcher searcher;
};


struct GlobalConfig
{
    BGIQD::stLFR::GraphEA_withBarcode graph_eab;
    BGIQD::stLFR::ContigCluster clusters;
    BGIQD::stLFR::ContigRoads roads;

    BGIQD::stLFR::ContigRoads extras;

    BGIQD::LOG::logger lger;
    BGIQD::SOAP2::FileNames fNames;
    enum FillStrategy
    {
        Unknow = 0 ,
        ShortestPath = 1,
        BarcodeCov_GiveUpCircle = 2 ,
        BarcodeCov_CutCircle = 3,
        BarcodeCov_FillCircle = 4,
    };
    FillStrategy strategy;
    float Ecov ;
    int max_length;
    int K;
    void Init(const std::string & prefix)
    {
        fNames.Init(prefix) ;
    }

    int thread;

    void FillContigRoad( int i ) //BGIQD::stLFR::ContigRoad & road)
    {
        auto & road = roads.roads[i];
        FillContigRoad1(road);
    }

    std::map<unsigned int , BGIQD::SOAP2::KeyEdge> edges;

    void LoadShortestPath()
    {
        std::string line ;
        int id  = 0 ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.connInfo());
        if(in == NULL )
            FATAL("failed to open xxx.connInfo for read !!! ");
        while( ! std::getline( *in , line).eof() )
        {
            auto items1 = BGIQD::STRING::split( line , "\t") ;
            assert( items1.size() >= 3 );
            unsigned int key = std::stoul(items1[0]);

            if( edges.find( key ) == edges.end() )
            {
                edges[key].Init( ++id , key );
            }

            if( items1[1] == "+" )
            {
                for(size_t i = 2 ; i < items1.size() ; i++ )
                {
                    edges.at(key).InitTo(items1[i]);
                }
            }
            else
            {
                for(size_t i = 2 ; i < items1.size() ; i++ )
                {
                    edges.at(key).InitFrom(items1[i]);
                }
            }
        }
        delete in ;
    }

    void ConstructRoadFill()
    {
        for( auto & a_road : roads.roads )
        {
            if( ! a_road.needMerge() )
                continue ;

            a_road.status = BGIQD::stLFR::ContigRoad::FillStatus::None ;
            a_road.fill_num = 0;
            for(int i=1  ; i < a_road.linear_length ; i++ )
            {
                auto start = a_road.getLinearStep(i);
                if( i == 1 )
                    a_road.contig_path.push_back(start.first) ;
                if( edges.find( start.first ) == edges.end() )
                {  // A' -->  X . but log as  X'-->A in from map of A
                    assert( edges.find( start.first -1 ) != edges.end() );

                    const auto & start_edge = edges.at(start.first - 1);
                    unsigned int  X1 = graph_eab.graph_ea.edge_array[start.second].bal_id ;

                    if( start_edge.from.find( X1 ) != start_edge.from.end() )
                    {
                        a_road.contig_path.insert( a_road.contig_path.end()
                                , start_edge.from.at(X1).path.begin()
                                , start_edge.from.at(X1).path.end() );
                    }
                    else
                    {
                        assert(0);
                    }
                }
                else 
                { // A -- > X . log as A -- > X in to map of A 
                    const auto & start_edge = edges.at(start.first);
                    if ( start_edge.to.find( start.second ) 
                            == start_edge.to.end() )
                    {
                        assert(0);
                    }
                    if ( ! start_edge.to.at(start.second).path.empty() )
                        a_road.contig_path.insert( a_road.contig_path.end()
                                , start_edge.to.at(start.second).path.begin()
                                , start_edge.to.at(start.second).path.end() );
                }

                a_road.contig_path.push_back(start.second) ;
                a_road.fill_num ++ ;
            }
            a_road.status = BGIQD::stLFR::ContigRoad::FillStatus::Complete;
        }
    }

    void RunAllJobs()
    {
        if( strategy == GlobalConfig::FillStrategy::ShortestPath )
        {
            //LoadShortestPath();
            ConstructRoadFill();
        }
        else
        {
            BGIQD::MultiThread::MultiThread t_jobs;
            t_jobs.Start(thread);

            for( int i= 0 ; i<(int)roads.roads.size(); i++ )
            {
                t_jobs.AddJob([this ,i](){ FillContigRoad(i);});
            }
            t_jobs.End();
            t_jobs.WaitingStop();
        }
        lger<<BGIQD::LOG::lstart()<<"fill contig road end ... "<<BGIQD::LOG::lend();
    }

    void FillContigRoad1( BGIQD::stLFR::ContigRoad & road )
    {
        static std::mutex write_mutex;
        static std::mutex extra_mutex;
        static std::mutex path_num_mutex;
        if ( ! road.needMerge() )
            return ;

        road.status = BGIQD::stLFR::ContigRoad::FillStatus::None ;
        road.fill_num = 0;
        if( strategy == GlobalConfig::FillStrategy::Unknow )
        {
            assert(0);
        }
        for(int i=1  ; i < road.linear_length ; i++ )
        {
            auto start = road.getLinearStep(i);
            SearchResult ret ;
            SearchAllPath(start.first , start.second, ret);
            if( ! ret.succ )
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
                if( road.fill_num == 0)
                {
                    road.contig_path.clear();
                    road.contig_path.push_back(ret.true_from);
                }
                if( strategy != GlobalConfig::FillStrategy::ShortestPath )
                {
                    BGIQD::stLFR::P2PGraph p2pgrapg;

                    p2pgrapg.base_graph = &graph_eab;
                    if( strategy  == GlobalConfig::FillStrategy::BarcodeCov_GiveUpCircle )
                        p2pgrapg.deal_circle = BGIQD::stLFR::P2PGraph::CircleStrategy::GiveUp ;
                    else if (strategy  == GlobalConfig::FillStrategy::BarcodeCov_CutCircle )
                        p2pgrapg.deal_circle = BGIQD::stLFR::P2PGraph::CircleStrategy::IgnoreCircle;
                    else if (strategy  == GlobalConfig::FillStrategy::BarcodeCov_FillCircle)
                        p2pgrapg.deal_circle = BGIQD::stLFR::P2PGraph::CircleStrategy::FillCircle ;
                    else 
                        assert(0);

                    p2pgrapg.ecov = Ecov ;
                    p2pgrapg.K =K ;

                    FindCorrectPath(ret.true_from , ret.true_to , ret , p2pgrapg );
                    // check if allpath find a correct path ?
                    {
                        std::lock_guard<std::mutex> l(path_num_mutex);
                        path_num_freq.Touch( p2pgrapg.path_num);
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
                    if( ! AppendPath( p2pgrapg ,ret,road ) )
                    {
                        road.status = BGIQD::stLFR::ContigRoad::FillStatus::Conflict ;
                        break;
                    }
                }
                else
                {
                    assert(0);
                }
                /*
                else
                {
                    auto extractShortestPath = [] (const Searcher & searcher , unsigned int from , unsigned int to)
                    {
                        unsigned int curr = to ;
                        std::stack<unsigned int> path;
                        while( curr != from )
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
                    auto retp = extractShortestPath( ret.searcher , ret.true_from , ret.true_to );
                    road.contig_path.insert(road.contig_path.end() 
                            , retp.begin()
                            , retp.end() );
                }*/
                road.status = BGIQD::stLFR::ContigRoad::FillStatus::PartSucc ;
                road.fill_num ++ ;
            }
        }

        if ( road.fill_num == road.linear_length - 1 )
        {
            road.status = BGIQD::stLFR::ContigRoad::FillStatus::Complete;
        }
        else
        {
            if( road.fill_num > 0 && road.fill_num < road.linear_length - 2 )
            {
                auto extra= road.Left(road.fill_num+1);
                FillContigRoad1(extra);
                {
                    std::lock_guard<std::mutex> l(extra_mutex);
                    extras.roads.push_back(extra);
                }
            }
        }
        {
            std::lock_guard<std::mutex> l(write_mutex);
            if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Conflict)
            {
                road_fill_freq.Touch("Conflict");
                return;
            }
            if( road.status == BGIQD::stLFR::ContigRoad::FillStatus::None )
            {
                road_fill_freq.Touch("None");
                return;
            }
            if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Complete)
            {
                road_fill_freq.Touch("Complete");
            }
            if ( road.status == BGIQD::stLFR::ContigRoad::FillStatus::PartSucc)
            {
                road_fill_freq.Touch("PartSucc");
            }
        }
    }
    //std::string updateEdge;
    //std::string arc;
    //std::string cluster;
    //std::string road;
    //std::string read2contig;
    //std::string contigroadfill;
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

    void LoadDatas()
    {
        roads.LoadRoads(fNames.contigroad());
        lger<<BGIQD::LOG::lstart()<<"load road end ... "<<BGIQD::LOG::lend();

        if( strategy == GlobalConfig::FillStrategy::ShortestPath )
        {
            graph_eab.graph_ea.LoadEdge(fNames.updatedEdge(),K);
            LoadShortestPath();
            lger<<BGIQD::LOG::lstart()<<"load connInfo end ... "<<BGIQD::LOG::lend();
        }
        else
        {
            //step1 .Loading files...
            graph_eab.graph_ea.LoadEdge(fNames.updatedEdge(),K);
            lger<<BGIQD::LOG::lstart()<<"load updateEdg end ... "<<BGIQD::LOG::lend();
            graph_eab.graph_ea.LoadArc(fNames.Arc());
            lger<<BGIQD::LOG::lstart()<<"load arc  end ... "<<BGIQD::LOG::lend();
            graph_eab.LoadBarcodeOnConfig(fNames.BarcodeOnContig());
            lger<<BGIQD::LOG::lstart()<<"load read2contig end ... "<<BGIQD::LOG::lend();
            clusters.loadCluster(fNames.cluster());
            lger<<BGIQD::LOG::lstart()<<"load cluster end ... "<<BGIQD::LOG::lend();
            LoadKeyInfo();
            lger<<BGIQD::LOG::lstart()<<"load key end ... "<<BGIQD::LOG::lend();
        }
    }
    void PrintRoadFill()
    {
        int filled = 0 ;
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.contigroadfill());
        if( out == NULL )
            FATAL("failed to open xxx.contigroadfill for write ");
        auto report_road = [&]( const BGIQD::stLFR::ContigRoad & road )
        {
            if( ! road.needMerge() )
                return;
            if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Conflict)
            {
                road_fill_freq.Touch("Conflict");
                return;
            }
            if( road.status == BGIQD::stLFR::ContigRoad::FillStatus::None )
            {
                road_fill_freq.Touch("None");
                return;
            }
            if (road.status == BGIQD::stLFR::ContigRoad::FillStatus::Complete)
            {
                road_fill_freq.Touch("Complete");
            }
            if ( road.status == BGIQD::stLFR::ContigRoad::FillStatus::PartSucc)
            {
                road_fill_freq.Touch("PartSucc");
            }
            bool has_circle = false ;
            for( auto i : road.circle_runs )
            {
                circle_run.Touch(i);
                if( i > 0 )
                    has_circle = true ;
            }
            if( has_circle )
            {
                *out<<"c\t";
            }
            else
            {
                *out<<"l\t";
            }
            for( const auto i : road.contig_path )
            {
                *out<<i<<'\t';
            }
            filled += (road.fill_num + 1 );
            *out<<std::endl;
            return ;
        };

        for( const auto &road : roads.roads)
        {
            report_road(road);
        }
        for( const auto & road: extras.roads )
        {
            report_road(road);
        }
        delete out ;
        lger<<BGIQD::LOG::lstart()<<"road fill use seed contig \n"<<filled<<BGIQD::LOG::lend();
    }
    

    void PrintStatistics()
    {
        lger<<BGIQD::LOG::lstart()<<"all path freq \n"<<path_num_freq.ToString()<<BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"road fill freq \n"<<road_fill_freq.ToString()<<BGIQD::LOG::lend();
        lger<<BGIQD::LOG::lstart()<<"circled  freq \n"<<circle_run.ToString()<<BGIQD::LOG::lend();
    }

    void SearchAllPath(unsigned int from  , unsigned int to , SearchResult & ret){

        // fill basic varibales
        BGIQD::SOAP2::Edge * edge_array = graph_eab.graph_ea.edge_array;
        ret.true_from = from ;
        ret.true_to = to ;

        std::map<unsigned int ,float> neibs ;
        auto & node = edge_array[to];
        unsigned int to_key_id ;
        if( node.IsKey() )
            to_key_id = to ;
        else
            to_key_id = node.bal_id ;
        neibs[to_key_id] = 0.5f;

        // Search graph
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

        ret.searcher.accesser.base = &graph_eab.graph_ea;
        ret.searcher.accesser.K = K ;
        ret.searcher.ender.Init( key , max_length,max_branch);
        ret.searcher.DoSPFSearch(ret.true_from);

        // Check result
        auto & ender = ret.searcher.ender ;
        ret.succ = true ;
        auto itr = ender.founder.find(to_key_id);
        if( itr == ender.founder.end() )
        {
            //assert(0);
            ret.succ = false;
        }
        auto & tos = itr->second;
        if( node.IsKey() &&  ! tos.base )
        {
            //assert(0);
            ret.succ = false;
        }
        else if ( ! node.IsKey() && ! tos.bal )
        {
            //assert(0);
            ret.succ = false;
        }
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
                    p2pgrapg.AddFromTo( a_to.prev , a_to.Base.value );
                    if( history.find( a_to.prev ) == history.end() )
                    {
                        history.insert( a_to.prev ) ;
                        nexts.push( result.searcher.fib_nodes.at(a_to.prev));
                    }
                }
                for( auto a_from : a_to.other_from )
                {
                    p2pgrapg.AddFromTo(a_from , a_to.Base.value ) ;
                    if( history.find( a_from ) == history.end() )
                    {
                        history.insert( a_from ) ;
                        nexts.push( result.searcher.fib_nodes.at(a_from));
                    }
                }
                nexts.pop();
            }
        }
        catch( ... )
        {
            assert(0);
        }

        p2pgrapg.RemovePalindromeLink();
        p2pgrapg.GeneratePath();
    }
    bool AppendPath( const  BGIQD::stLFR::P2PGraph & p2pgrapg , const SearchResult & result, BGIQD::stLFR::ContigRoad & road)
    {
        if( ! road.needMerge() )
            return false;
        if( ! result.succ )
            return false;
        if( p2pgrapg.path_num < 1 )
            return false;

        road.contig_path.insert(road.contig_path.end() 
                , p2pgrapg.final_path.begin()
                , p2pgrapg.final_path.end() );
        road.contig_path.push_back(result.true_to);

        road.circle_runs.push_back(p2pgrapg.final_circled);
        return true;
    }
    int max_branch;
} config;





int  main(int argc, char **argv)
{
    //step0 Parse parmeters...
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix \n\
                                                  Input xxx.Arc\n\
                                                        xxx.updated.edge\n\
                                                        xxx.cluster\n\
                                                        xxx.contigroad\n\
                                                        xxx.barcodeOnContig;\n\
                                                 Output xxx.contigroadfill");
    DEFINE_ARG_REQUIRED(int , kvalue, "K value");
    DEFINE_ARG_OPTIONAL(int , thread,"thread num ","8");
    DEFINE_ARG_OPTIONAL(float , Ecov, "Ecov of contigs. must set this if fill circle.","10.0");
    DEFINE_ARG_REQUIRED(int, fill_strategy, "fill strategy \n\
                                                        ShortestPath = 1\n\
                                                        BarcodeCov_GiveUpCircle = 2\n\
                                                        BarcodeCov_CutCircle = 3\n\
                                                        BarcodeCov_FillCircle = 4\
    ");
    DEFINE_ARG_OPTIONAL(int, searchDepth, "search depth (bp) ","10000");
    DEFINE_ARG_OPTIONAL(int, maxBranch,"max search branch ","10");
    END_PARSE_ARGS

    BGIQD::LOG::timer t(config.lger,"FillContigRoad");
    config.max_branch = maxBranch.to_int();
    config.K = kvalue.to_int();
    config.Init(prefix.to_string());
    config.max_length = searchDepth.to_int();
    config.Ecov = Ecov.to_float();
    config.strategy = static_cast<GlobalConfig::FillStrategy>(fill_strategy.to_int());
    config.thread = thread.to_int();
    BGIQD::LOG::logfilter::singleton().get("FillContigRoad",BGIQD::LOG::loglevel::INFO , config.lger);
    config.lger<<BGIQD::LOG::lstart()<<"parse args end ... "<<BGIQD::LOG::lend();

    //step1 Load data from disk...
    config.LoadDatas() ;
    //step2 Fill roads ...
    config.RunAllJobs();
    //step3 Print ...
    config.PrintRoadFill();
    //step4 Print log
    config.PrintStatistics();
}
