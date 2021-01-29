#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/string/stringtools.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"

#include "utils/misc/fileName.h"
#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include "utils/linear_fitting/Minimum_multiplication.h"

#include <map>
#include <vector>
#include <array>
#include <algorithm>

struct AppConfig
{
    BGIQD::LOG::logger loger;
    BGIQD::MISC::FileNames fName;

    typedef BGIQD::LINEARFITTING::Item<int,float> item;
    typedef BGIQD::LINEARFITTING::Linear<int,float> LC;
    typedef BGIQD::stLFR::TrunkGap<int> GapInfo;

    std::map<int,std::vector<GapInfo>> gaps;
    BGIQD::stLFR::BinRelationArray  bra;
    std::map<unsigned int , float > to_next_sim ;

    std::vector<item> gap_real_data;

    LC lc;

    //  ( contig_id , binId ) <--> ( start_pos , end_pos )
    std::map<std::pair<unsigned int , int > , std::pair<int,int> > bin_on_contig_pos ;

    void CacheBinPos()
    {
        for(int i = 0 ; i <(int) bra.size() ; i++ )
        {
            const auto & bin = bra[i];
            unsigned int contig = bin.contigId ;
            int binId = bin.binId ;
            bin_on_contig_pos[std::make_pair( contig , binId) ] 
                = std::make_pair( bin.start ,bin.end ) ;
        }
    }

    bool IsA2BNext(unsigned int A , unsigned int B )
    {
        if( contig_indexs.find( A ) == contig_indexs.end()) 
            return false ;
        if( contig_indexs.find( B ) == contig_indexs.end()) 
            return false ;
        return contig_indexs.at(B) - contig_indexs.at(A) == 1 ;
    }

    void LoadNextSim()
    {
        for(int i = 0 ; i <(int) bra.size() ; i++ )
        {
            const auto & bin = bra[i];
            unsigned int contig = bin.contigId ;
            for( const auto & pair : bin.sims)
            {
                const auto & next = pair.second;
                if(!  IsA2BNext( contig , next.contigId ) )
                    continue ;
                if( to_next_sim.find( contig ) ==to_next_sim.end() )
                    to_next_sim[contig]  = next.simularity ;
                else if(to_next_sim[contig] < next.simularity)
                    to_next_sim[contig] = next.simularity ;
            }
        }
    }
    void LoadLinearCachce()
    {
        for(int i = 0 ; i <(int) bra.size() ; i++ )
        {
            const auto & bin = bra[i];
            unsigned int contig = bin.contigId ;
            int binId = bin.binId ;
            for( const auto & pair : bin.sims)
            {
                const auto & next = pair.second;
                if( next.contigId == contig && next.binId > binId )
                {
                    const auto & this_pos = bin_on_contig_pos.at( 
                            std::make_pair( contig , binId) ) ;
                    const auto & other_pos = bin_on_contig_pos.at( 
                            std::make_pair( contig , next.binId) ) ;
                    int gap = 0 ;
                    if( this_pos.first < other_pos.first )
                        gap = other_pos.first - this_pos.second ;
                    else
                        gap = this_pos.first - other_pos.second ;

                    gap_real_data.push_back( item { gap , next.simularity} );
                }
            }
        }
        lc = BGIQD::LINEARFITTING::lineFit(gap_real_data);
    }

    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
        bra.Init(10000);
    }

    void LoadTrunk()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        BGIQD::stLFR::Load_MST_Trunk_Linear(*in, gaps);
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load Trunk done "<<BGIQD::LOG::lend() ;
    }

    void  LoadBinRelationArrayFromFile()
    {
        BGIQD::stLFR::LoadBinRelationArray(fName.bin_cluster("gap"),bra);
        loger<<BGIQD::LOG::lstart() << "Load Bin cluster done "<<BGIQD::LOG::lend() ;
    }

    // contig_id <--> index , used to find neighbors.
    std::map<int , int > contig_indexs ;

    void BuildNeighbors()
    {
        int i = 0 ;
        for( const auto & pair : gaps )
        {
            i ++ ;
            const auto & a_scaff = pair.second ;
            bool start = true ;
            for( const auto & a_gap : a_scaff )
            {
                if( start )
                {
                    i++ ;
                    contig_indexs[a_gap.prev] = i ; 
                    start = false ;
                }
                i ++ ;
                contig_indexs[a_gap.next] = i ;
            }
        }
    }

    void PrintGapSim()
    {
        auto out = BGIQD::FILES::FileWriterFactory::
            GenerateWriterFromFileName(fName.gap_sim());
        if( out == NULL )
            FATAL( "failed to open xxx.gap_sim to write ");
        for( const auto & pair :  to_next_sim )
        {
            (*out)<<pair.first<<'\t'
                <<int(pair.second*1000000)
                <<'\n';
        }
        delete out ;
    }

    void PrintGapArea()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.gap_area());
        if( out == NULL )
            FATAL( "failed to open xxx.gap_area to write ");
        for( int i = 1 ; i < max_gap ; i++ )
        {
            (*out)<<i<<'\t'<<lc.getY(i)<<'\n';
        }
        delete out ;
        return ;
    }
    int max_gap ;

} config;

int main(int argc, char **argv)
{
    //step 0 Parse parmeters...
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix of read names.\n\
                                                    In \n\
                                                        xxx.mintree_trunk_linear; xxx.bin_cluster ;\n\
                                                    Out\n\
                                                        xxx.gap_sim ; xxx.gap_area");
    DEFINE_ARG_OPTIONAL(int , max_gap , "the max gap size ", "15000"); 
    END_PARSE_ARGS;

    config.max_gap = max_gap.to_int() ;
    config.Init( prefix.to_string());
    config.LoadTrunk();
    config.BuildNeighbors();
    config.LoadBinRelationArrayFromFile();
    config.LoadNextSim();
    config.PrintGapSim();
    config.CacheBinPos() ;
    config.LoadLinearCachce();
    config.PrintGapArea();
    return 0;
}
