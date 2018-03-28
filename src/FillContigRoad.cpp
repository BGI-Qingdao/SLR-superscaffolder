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
    BGIQD::FREQ::Freq<std::string> road_fill_freq;
} config;

struct DepthSearchResult
{
    std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > paths;
    std::map<unsigned int , std::vector<std::list<BGIQD::SOAP2::Edge> > > mids;
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
};

DepthSearchResult SearchAllPath(unsigned int from  , unsigned int to){

    DepthSearchResult ret;
    ret.head_tail = true ;
    ret.downstream = true ;
    ret.status = DepthSearchResult::NotRearch;
    BGIQD::SOAP2::Edge * edge_array = config.graph_eab.graph_ea.edge_array;

    unsigned int head_id = from;
    unsigned int tail_id = to;
    unsigned int search_id = head_id;
    unsigned int target_id = tail_id;

    /*
    if( edge_array[head_id].bal_id == head_id && edge_array[tail_id].bal_id != tail_id)
    {
        ret.head_tail = !ret.head_tail ;
        unsigned int tmp = search_id ;
        search_id = tail_id;
        target_id = tmp;
    }
    */

    std::list<BGIQD::SOAP2::Edge> stack;
    std::map<unsigned int , BGIQD::SOAP2::Edge > history;
    std::map<unsigned int ,float> neibs ;
    neibs[target_id] = 1.0f;
    // try downstream 
    edge_array[search_id].DepthSearch(edge_array,stack,
            history, ret.paths , ret.mids ,edge_array[search_id].length , neibs);
    /*
    if( ret.paths.empty() )
    {
        if( edge_array[search_id].bal_id  == search_id )
        {
            // now both contig has no bal_id.
            // try another downstream
            ret.head_tail = !ret.head_tail ;
            unsigned int tmp = search_id ;
            search_id = tail_id;
            target_id = tmp;
            neibs.clear();
            neibs[target_id] = 1.0f;
        }
        else 
        {
            //try upstream
            search_id = edge_array[search_id].bal_id;
            ret.downstream = false ;
        }
        stack.clear();
        history.clear();
        ret.paths.clear();
        ret.mids.clear();

        edge_array[search_id].DepthSearch(edge_array,stack,
                history, ret.paths , ret.mids ,edge_array[search_id].length , neibs);
    }
    */

    if( ret.paths.empty() )
    {
        ret.status = DepthSearchResult::NotRearch;
        return ret;
    }
    //detect relationship
    unsigned int true_target = (ret.paths.begin()->second)[0].front().id;

    ret.true_from = search_id ;
    ret.true_to = true_target ;
    assert( true_target == target_id );
    ret.status = DepthSearchResult::A1B1_B2A2 ;

    /*
    if(  ret.head_tail ) 
    {
        // Start from A contig
        if( ret.downstream )
            if( true_target == target_id )
                ret.status = DepthSearchResult::A1B1_B2A2 ;
            else 
                ret.status = DepthSearchResult::A1B2_B1A2 ;
        else 
            if( true_target == target_id )
                ret.status = DepthSearchResult::B2A1_A2B1 ;
            else 
                ret.status = DepthSearchResult::B1A1_A2B2 ;
    }
    else
    {
        // Start from B contig
        if( ret.downstream )
            if( true_target == target_id )
                ret.status = DepthSearchResult::B1A1_A2B2 ;
            else 
                ret.status = DepthSearchResult::A1B2_B1A2 ;
        else 
            if( true_target == target_id )
                ret.status = DepthSearchResult::B2A1_A2B1 ;
            else 
                ret.status = DepthSearchResult::A1B1_B2A2 ;
    }
    */
    return ret ;
}

void FindCorrectPath(unsigned int from , unsigned int to, const DepthSearchResult & result , BGIQD::stLFR::P2PGraph & p2pgrapg )
{
    p2pgrapg.Init(from,to);
    for( const auto & i : result.paths )
        p2pgrapg.AddPath(to , i.second);
    for( const auto & i : result.mids )
        p2pgrapg.AddMid(i.first , i.second);
    p2pgrapg.GeneratePath();
}

bool AppendPath( const  BGIQD::stLFR::P2PGraph & p2pgrapg , const DepthSearchResult & result, BGIQD::stLFR::ContigRoad & road)
{
    if(! road.needMerge() )
        return false;
    if( result.status== DepthSearchResult::NotRearch )
        return false;
    if( p2pgrapg.path_num < 1 )
        return false;

    if( result.head_tail )
    {
        road.contig_path.insert(road.contig_path.end() 
                , p2pgrapg.final_path.paths.begin()
                , p2pgrapg.final_path.paths.end() );
    }
    else
    {
        road.contig_path.insert(road.contig_path.end() 
                , p2pgrapg.final_path.paths.rbegin()
                , p2pgrapg.final_path.paths.rend() );
    }

    return true;
}

void FillContigRoad( BGIQD::stLFR::ContigRoad & road)
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
        auto ret = SearchAllPath(start.first , start.second);
        if( ret.status == DepthSearchResult::NotRearch )
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
            if ( ret.status == DepthSearchResult::A1B1_B2A2 ||  ret.status == DepthSearchResult::A1B2_B1A2  )
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
            BGIQD::stLFR::P2PGraph p2pgrapg;
            p2pgrapg.base_graph = &config.graph_eab;
            FindCorrectPath(ret.true_from , ret.true_to , ret , p2pgrapg);
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
            if( road.fill_num== 0 )
            {
                if( ret.head_tail )
                    road.contig_path.push_back( ret.true_from );
                else 
                    road.contig_path.push_back( ret.true_to );
            }
            if( ! AppendPath( p2pgrapg ,ret,road ) )
            {
                road.status = BGIQD::stLFR::ContigRoad::FillStatus::Conflict ;
                break;
            }
            if( ret.head_tail )
                road.contig_path.push_back( ret.true_to );
            else 
                road.contig_path.push_back( ret.true_from );
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
    END_PARSE_ARGS

    config.K = kvalue.to_int();
    config.arc = prefix.to_string() +".Arc";
    config.updateEdge = prefix.to_string() +".updated.edge";
    config.cluster= prefix.to_string() +".cluster";
    config.road = prefix.to_string() +".contigroad";
    config.read2contig = prefix.to_string() +".read2contig";

    if(! t_num.setted )
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

        for(int i= 0 ; i<(int)config.roads.roads.size(); i++)
        {
            t_jobs.AddJob([i](){ FillContigRoad(std::ref(config.roads.roads[i])); });
        }
        t_jobs.End();
        t_jobs.WaitingStop();
    }

    config.lger<<BGIQD::LOG::lstart()<<"fill contig road end ... "<<BGIQD::LOG::lend();

    report();

    config.lger<<BGIQD::LOG::lstart()<<"report end ... "<<BGIQD::LOG::lend();
    config.lger<<BGIQD::LOG::lstart()<<"all path freq \n"<<config.path_num_freq.ToString()<<BGIQD::LOG::lend();
    config.lger<<BGIQD::LOG::lstart()<<"road fill freq \n"<<config.road_fill_freq.ToString()<<BGIQD::LOG::lend();

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
