#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "algorithm/graph/Graph.h"
#include "algorithm/graph/DepthSearch.h"

#include "soap2/fileName.h"
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
        std::set<unsigned int> relations;
        void InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            ist>>prev>>next;
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
    BGIQD::SOAP2::FileNames fName;
    int searchMax ;
    int total_fill ; 
    std::set<unsigned int> trunk_seeds;
    std::vector<GapInfo> all_gaps;
    BGIQD::stLFR::ContigPEGraph pe_graph;
    std::map<unsigned int , int > contigLen_cache;

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
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds_cluster_seeds());
        if( in == NULL )
            FATAL( " failed to open xxx.seeds_cluster_seeds to read !!!" );

        auto parseline = [this](const std::string & line) -> void 
        {
            if( line == "bigraph {" || line == "}" )
                return ;
            BGIQD::stLFR::PEEdge edge;
            edge.InitFromString(line);
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

        for (auto & gap : all_gaps)
        {
            BGIQD::stLFR::ContigPEGraph sub_graph = pe_graph.SubGraph(gap.relations);
            for(auto & node : sub_graph.nodes )
            {
                sub_graph.MakeEdgeNext(node.second.id);
            }
            int succ = 0;
            //try 1 order
            {
                typedef BGIQD::GRAPH::EdgeIterator<BGIQD::stLFR::ContigPEGraphAccess> EdgeItr;

                typedef BGIQD::GRAPH::DepthSearch<
                    BGIQD::stLFR::ContigPEGraphAccess,
                    EdgeItr,
                    BGIQD::stLFR::DepthEnder
                        > Searcher;
                Searcher searcher;
                searcher.accesser.base = &sub_graph;
                searcher.ender.Init(typer, searchMax );
                searcher.DoDepthSearch(gap.prev,0);
                if( searcher.nodes.find(gap.next) != searcher.nodes.end() 
                        || searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    succ ++ ;
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
                searcher.ender.Init(typer, searchMax );
                if( searcher.nodes.find(gap.next) != searcher.nodes.end() 
                        || searcher.nodes.find(gap.next +1 ) != searcher.nodes.end() )
                {
                    succ ++ ;
                }
            }

            if ( succ > 0 )
                total_fill ++ ;

        }
    }


    void PrintFillResult()
    {
        std::cout<<"total fill "<<total_fill<<" in "<<all_gaps.size()<<std::endl;
    }

    void Init(const std::string & prefix , int search_len_max )
    {
        fName.Init(prefix);
        searchMax = search_len_max ;
        BGIQD::LOG::logfilter::singleton().get("FillTrunkByPE", BGIQD::LOG::loglevel::INFO,loger);
        total_fill = 0 ;
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.seeds()) ;
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
    END_PARSE_ARGS;
    config.Init(prefix.to_string() , searchMax.to_int());
    BGIQD::LOG::timer t(config.loger,"FillTrunkByPE");
    config.LoadSeedCluster();
    config.LoadSeeds();
    config.LoadPEGraph();
    config.CalcAll();
    config.PrintFillResult();

    return 0;
}
