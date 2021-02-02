#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/string/stringtools.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"

#include "utils/graph/DepthFirstSearch.h"

#include "utils/misc/fileName.h"
#include "stLFR/contigPEGraph.h"
#include "stLFR/contigPEGraphSearch.h"

#include <set>
#include <sstream>

struct AppConfig
{
    struct GapInfo
    {
        unsigned int prev ;
        unsigned int next ;
        unsigned int prev_1 ;
        unsigned int next_1 ;
        std::set<unsigned int> relations;
        void InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            ist>>prev>>next;
            relations.insert(prev);
            relations.insert(prev+1);
            relations.insert(next);
            relations.insert(next+1);
            while( !ist.eof() )
            {
                unsigned int r ;
                ist>>r ;
                relations.insert(r);
                relations.insert(r+1);
            }
        }
        std::vector<unsigned int> path;
    };

    BGIQD::LOG::logger loger;
    BGIQD::MISC::FileNames fName;
    int searchMax ;
    int total_fill ; 
    std::set<unsigned int> trunk_seeds;
    std::vector<GapInfo> all_gaps;
    BGIQD::stLFR::ContigPEGraph pe_graph;
    std::map<unsigned int , int > contigLen_cache;
    BGIQD::MISC::Freq<int> fill_freq;
    int is_max;
    int min_count ;
    bool use_all_pe_graph;
    void LoadSeedCluster()
    {
        auto parseline = [this](const std::string & line) -> void 
        {
            GapInfo tmp;
            tmp.InitFromString(line);
            trunk_seeds.insert(tmp.prev);
            trunk_seeds.insert(tmp.next);
            all_gaps.push_back(tmp);
        };
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds_cluster_seeds());
        if( in == NULL )
            FATAL( " failed to open xxx.seeds_cluster_seeds to read !!!" );

        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
    }

    void LoadPEGraph()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.pe_graph());
        if( in == NULL )
            FATAL( " failed to open xxx.pe_graph to read !!!" );

        auto parseline = [this](const std::string & line) -> void 
        {
            if( line == "digraph {" || line == "}" )
                return ;
            BGIQD::stLFR::PEEdge edge;
            edge.InitFromString(line);
            if( std::abs(edge.len) > is_max )
                return ;
            if( edge.count < min_count ) 
                return ;
            pe_graph.AddEdge(edge);
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete  in;
    }

    void CalcAll()
    {

        auto typer = [this](unsigned int id )
        {
            if( trunk_seeds.find(id) != trunk_seeds.end() )
                return BGIQD::stLFR::DepthEnder::NodeType::Seed ;
            else
                return BGIQD::stLFR::DepthEnder::NodeType::Others;
        };

        typedef BGIQD::GRAPH::EdgeIterator<BGIQD::stLFR::ContigPEGraphAccess> EdgeItr;

        typedef BGIQD::GRAPH::DepthSearch<
            BGIQD::stLFR::ContigPEGraphAccess,
            EdgeItr,
            BGIQD::stLFR::DepthEnder
                > Searcher;

        auto get_path = [](const Searcher & searcher , unsigned int root , unsigned int end )
        {
            std::vector<unsigned int> ret ;
            if( searcher.nodes.find( end) == searcher.nodes.end() )
                return ret;
            std::stack<unsigned int> path;
            auto & ender = searcher.nodes.at( end ) ;
            unsigned int prev = ender.prev ;
            int max_depth = 1000 ;
            int curr = 0 ;
            path.push( end );
            while( prev != root  && curr < max_depth )
            {
                path.push(prev) ;
                auto & prever = searcher.nodes.at(prev) ;
                prev = prever.prev ;
                curr ++ ;
            }
            if( prev== root )
            {
                path.push(root);
                while(!path.empty() )
                {
                    ret.push_back( path.top() );
                    path.pop();
                }
            }
            return ret;
        };

        for(auto & node : pe_graph.nodes )
        {
            pe_graph.MakeEdgeNext(node.second.id);
        }

        for (auto & gap : all_gaps)
        {
            BGIQD::stLFR::ContigPEGraph cluster_sub_graph;
            BGIQD::stLFR::ContigPEGraph * sub_graph_ptr ;
            if( !use_all_pe_graph )
            {
                cluster_sub_graph = pe_graph.SubGraph(gap.relations);
                sub_graph_ptr = &cluster_sub_graph ;
            }
            else 
                sub_graph_ptr = &pe_graph ;
            BGIQD::stLFR::ContigPEGraph & sub_graph = *(sub_graph_ptr);
            if( ! use_all_pe_graph )
            {
                for(auto & node : sub_graph.nodes )
                {
                    sub_graph.MakeEdgeNext(node.second.id);
                }
            }
            std::vector<std::vector<unsigned int> > paths;
            unsigned int s1 , e1 ;
            //try 1 order
            {
                Searcher searcher;
                searcher.accesser.base = &sub_graph;
                searcher.accesser.Init();
                searcher.ender.Init(typer, searchMax );
                searcher.DoDepthSearch(gap.prev,0);
                if( searcher.nodes.find(gap.next) != searcher.nodes.end()  )
                {
                    s1= gap.prev ;
                    e1 = gap.next;
                    paths.push_back(get_path(searcher,gap.prev,gap.next));
                }
                if( searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    paths.push_back(get_path(searcher,gap.prev,gap.next+1));
                    s1= gap.prev ;
                    e1 = gap.next +1;
                }
            }
            // try another order
            {
                typedef BGIQD::GRAPH::EdgeIterator<BGIQD::stLFR::ContigPEGraphAccess> EdgeItr;

                typedef BGIQD::GRAPH::DepthSearch<
                    BGIQD::stLFR::ContigPEGraphAccess,
                    EdgeItr,
                    BGIQD::stLFR::DepthEnder
                        > Searcher;
                Searcher searcher;
                searcher.accesser.base = &sub_graph;
                searcher.accesser.Init();
                searcher.ender.Init(typer, searchMax );
                searcher.DoDepthSearch(gap.prev+1,0);
                if( searcher.nodes.find(gap.next) != searcher.nodes.end()  )
                {
                    paths.push_back(get_path(searcher,gap.prev+1,gap.next));
                    s1= gap.prev +1 ;
                    e1 = gap.next;
                }
                if( searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    paths.push_back(get_path(searcher,gap.prev+1,gap.next+1));
                    s1= gap.prev + 1 ;
                    e1 = gap.next +1;
                }
            }
            fill_freq.Touch(paths.size());
            if( paths.size() == 1 )
            {
                gap.path = paths[0] ;
                gap.prev_1 = s1 ;
                gap.next_1 = e1 ;
            }
        }
        loger<<BGIQD::LOG::lstart()<<fill_freq.ToString()<<BGIQD::LOG::lend();
    }

    void PrintFillResult()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.trunk_fill());
        if( out == NULL )
            FATAL(" failed to open xxx.trunk_fill to write !! ");
        for( const auto & gap :all_gaps )
        {
            if( gap.path.empty() )
                continue ;

            (*out)<<gap.prev<<'\t'<<gap.next<<'\t'
                <<gap.prev_1<<'\t'<<gap.next_1;
            for( auto i : gap.path )
            {
                (*out)<<'\t'<<i;
            }
            (*out)<<'\n';
        }
        delete out ;
    }
    void PrintTest()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName("test.fill");
        if( out == NULL )
            FATAL(" failed to open test.fill to write !! ");
        for( const auto & gap :all_gaps )
        {
            if( gap.path.empty() )
                continue ;

            (*out)<<gap.path.size();
            for( auto i : gap.path )
            {
                (*out)<<'\t'<<i;
            }
            (*out)<<'\n';
        }
        delete out ;
    }

    void Init(const std::string & prefix , int search_len_max , int max_is , bool  use_pe)
    {
        fName.Init(prefix);
        searchMax = search_len_max ;
        BGIQD::LOG::logfilter::singleton().get("FillTrunkByPE", BGIQD::LOG::loglevel::INFO,loger);
        total_fill = 0 ;
        is_max = max_is ;
        use_all_pe_graph = use_pe;
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds("pe")) ;
        if( in == NULL )
            FATAL( "open .seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            auto items = BGIQD::STRING::split(line,"\t");
            contigLen_cache[std::stoul(items[0])] = std::stoul(items[1]);
            pe_graph.AddNode(std::stoul(items[0]) , std::stoul(items[1]));
        }
        delete in ;
    };

} config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string, prefix ,"prefix of files.");
        DEFINE_ARG_OPTIONAL(int, searchMax,"max search length." ,"5000");
        DEFINE_ARG_OPTIONAL(int, insert_max,"max insert_size avaliable." ,"1000");
        DEFINE_ARG_OPTIONAL(int, min_count ,"min valid count ." ,"2");
        DEFINE_ARG_OPTIONAL(bool,use_all_pe_graph,"use all pe_graph instead of cluster sub graph." ,"false");
        //DEFINE_ARG_OPTIONAL(bool,ptest,"print date for test." ,"false");
    END_PARSE_ARGS;
    config.Init(prefix.to_string() , searchMax.to_int(),insert_max.to_int(),use_all_pe_graph.to_bool());
    config.min_count = min_count.to_int();
    BGIQD::LOG::timer t(config.loger,"FillTrunkByPE");
    config.LoadSeedCluster();
    config.LoadSeeds();
    config.LoadPEGraph();
    config.CalcAll();
    config.PrintFillResult();
    //if( ptest.to_bool() )
    //    config.PrintTest();
    return 0;
}
