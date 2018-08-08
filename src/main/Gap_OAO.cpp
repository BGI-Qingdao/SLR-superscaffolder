#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"

#include "soap2/fileName.h"
#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <map>
#include <vector>
#include <array>
#include <algorithm>

struct AppConfig
{
    struct ScaffItem
    {
        unsigned int base ;
        std::vector<float> LLeft;
        std::vector<float> LRight;
        std::vector<float> RLeft;
        std::vector<float> RRight;

        unsigned int Lctg;
        int Lquality;
        unsigned int Rctg;
        int Rquality;
        void Calc()
        {
            Calc1Side( LLeft , LRight , Lctg , Lquality);
            Calc1Side( RLeft , RRight , Rctg , Rquality);
        }

        void Calc1Side( const std::vector<float> & left , const std::vector<float> & right , unsigned int & ctg , int & quality )
        {
            int l = 0,r = 0;
            float lt = 0, rt = 0;
            for( int i = 0 ; i < (int)left.size() ; i++)
            {
                lt += left[i] ;
                rt += right[i] ;
                if( left[i] < right[i] )
                    l ++ ;
                else
                    r++ ;
            }
            if( l < r ) 
            {
                ctg = base + 1;
                quality = l -r ;
            }
            else if ( r < l )
            {
                ctg = base ;
                quality = r -l ;
            }
            else
            {
                if( lt < rt )
                    ctg = base +1 ;
                else 
                    ctg = base ;
                quality = 0 ;
            }
        }
    };

    std::map<int,std::vector<ScaffItem>> scaffs;

    std::map<unsigned int , int[2]> contigIndex;
    int rank ;

    struct GapExtra
    {
        unsigned int true_prev ;
        unsigned int true_next ;
        std::array<float,4> bin_index ;
        float value;
        float sim ;
    };

    typedef BGIQD::stLFR::TrunkGap<GapExtra> GapInfo;
    std::map<int,std::vector<GapInfo>> gaps;
    BGIQD::stLFR::BinRelationArray  bra;
    BGIQD::LOG::logger loger;
    BGIQD::SOAP2::FileNames fName;

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
        BGIQD::stLFR::LoadBinRelationArray(fName.bin_cluster(),bra);
        loger<<BGIQD::LOG::lstart() << "Load Bin cluster done "<<BGIQD::LOG::lend() ;
    }

    void BuildSims()
    {
        // prepare contigIndex
        for(int i = 0 ; i <(int) bra.size() ; i++ )
        {
            const auto & bin = bra[i];
            unsigned int contig = bin.contigId;
            assert( bin.binId < 2 );
            contigIndex[contig][bin.binId] = i ;
        }

        // prepare scaffs
        for( auto & pairs: gaps )
        {
            for( auto & gap : pairs.second )
            {
                if ( scaffs[pairs.first].empty() )
                {
                    ScaffItem item ;
                    item.base = gap.prev ;
                    scaffs[pairs.first].push_back(item);
                }
                ScaffItem item ;
                item.base = gap.next;
                scaffs[pairs.first].push_back(item);
            }
        }
        // fill scaffs 
        for(auto & pair : scaffs)
        {
            for(size_t i = 0 ; i < pair.second.size() ; i++ )
            {
                int l = i - rank > 0 ?  i - rank : 0 ;
                int r = i + rank < pair.second.size() ? i + rank : pair.second.size() -1 ;
                auto & item = pair.second[i] ;
                unsigned int contig = item.base ;
                std::map<unsigned int , int > lindex;
                std::map<unsigned int , int > rindex;
                for ( size_t j = l ; j < i ; j ++ )
                {
                    lindex[pair.second[j].base]=item.LLeft.size();
                    item.LLeft.push_back(0);
                    item.LLeft.push_back(0);
                    item.LRight.push_back(0);
                    item.LRight.push_back(0);
                }
                for( size_t j = r ; j > i ; j -- )
                {
                    lindex[pair.second[j].base]=item.LLeft.size();
                    item.RLeft.push_back(0);
                    item.RLeft.push_back(0);
                    item.RRight.push_back(0);
                    item.RRight.push_back(0);
                }

                const auto & binLeft = bra[contigIndex[contig][0]];
                for( const auto & pair : binLeft.sims)
                {
                    const auto & sim = pair.second ;
                    if( lindex.find(sim.contigId) != lindex.end() )
                    {
                        int index = lindex[sim.contigId] + sim.binId;
                        item.LLeft[index] = sim.simularity ;
                    }
                    if( rindex.find(sim.contigId) != rindex.end() )
                    {
                        int index = rindex[sim.contigId] + sim.binId;
                        item.RLeft[index] = sim.simularity ;
                    }
                }

                const auto & binRight = bra[contigIndex[contig][1]];
                for( const auto & pair : binRight.sims)
                {
                    const auto & sim = pair.second ;
                    if( lindex.find(sim.contigId) != lindex.end() )
                    {
                        int index = lindex[sim.contigId] + sim.binId;
                        item.LRight[index] = sim.simularity ;
                    }
                    if( rindex.find(sim.contigId) != rindex.end() )
                    {
                        int index = rindex[sim.contigId] + sim.binId;
                        item.RRight[index] = sim.simularity ;
                    }
                }
            }
        }
    }
    void GetGapSim()
    {
        for( auto & pairs: gaps )
        {
            for( auto & gap : pairs.second )
            {
                gap.data.bin_index[0] = 0.000001;
                gap.data.bin_index[1] = 0.000001;
                gap.data.bin_index[2] = 0.000001;
                gap.data.bin_index[3] = 0.000001;
            }
        }

        for(int i = 0 ; i <(int) bra.size() ; i++ )
        {
            const auto & bin = bra[i];
            unsigned int contig = bin.contigId;
            assert( bin.binId < 2 );
            for( auto & pairs: gaps )
            {
                for( auto & gap : pairs.second )
                {
                    if( gap.prev == contig )
                    {
                        int start = bin.binId * 2 ;
                        for( const auto & pair : bin.sims)
                        {
                            const auto & sim = pair.second ;
                            if( sim.contigId == gap.next )
                            {
                                gap.data.bin_index[start+sim.binId] = sim.simularity;
                            }
                        }
                    }
                }
            }
        }
        loger<<BGIQD::LOG::lstart() << "GetGapSim done "<<BGIQD::LOG::lend() ;
    }
    void CalcAll1()
    {
        for( auto & pair: scaffs)
        {
            for( auto & item : pair.second)
            {
                item.Calc();
            }
        }

    }
    void CalcAll()
    {
        int total = 0 ;
        int failed = 0 ;
        for( auto & pairs: gaps )
        {
            for( auto & gap : pairs.second )
            {
                float max  = -0.1f;
                float min = 1.1f;
                for( float x :gap.data.bin_index) 
                {
                    if( x > max )
                        max = x ;
                    if( x < min )
                        min = x ;
                }
                assert( 0 < max && max  < 1 );
                assert( 0 < min && min < 1 );
                gap.data.value = max / min ;
                gap.data.sim = max ;
                /*
                 *   R1  L1   L2 R2
                 *   <-----   ----->
                 */
                // LL max , RR min
                if( max == gap.data.bin_index[0]) // && min == gap.data.bin_index[3] )
                {
                    gap.data.true_prev = gap.prev +1 ;
                    gap.data.true_next = gap.next ;
                }
                /*
                 *   L1  R1   R2 L2
                 *   ----->   <-----
                 */
                // RR  max , LL min 
                else if( max == gap.data.bin_index[3] )// && min == gap.data.bin_index[0] )
                {
                    gap.data.true_prev = gap.prev ;
                    gap.data.true_next = gap.next +1 ;
                }
                /*
                 *   L1  R1   L2 R2
                 *   ----->   ----->
                 */
                // RL max LR min
                else if( max == gap.data.bin_index[2] )//&& min == gap.data.bin_index[1] )
                {
                    gap.data.true_prev = gap.prev ;
                    gap.data.true_next = gap.next ;
                }
                /*
                 *   R1  L1   R2 L2
                 *   <-----   <-----
                 */
                // LR max RL min
                else if( max == gap.data.bin_index[1] )//&& min == gap.data.bin_index[2] )
                {
                    gap.data.true_prev = gap.prev +1 ;
                    gap.data.true_next = gap.next +1 ;
                }
                else
                {
                    gap.data.true_prev = 0;
                    gap.data.true_next = 0;
                    failed ++ ;
                }
                total ++ ;
            }
        }
        loger<<BGIQD::LOG::lstart() << "done , succ "<<total - failed << " in "<<total<<BGIQD::LOG::lend() ;
    }

    void PrintGapOO1()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.gap_oo());
        if( out == NULL )
            FATAL( " failed to open xxx.gap_oo for write ");

        for( auto & pair: scaffs)
        {
            ScaffItem * prev = NULL ;
            for( auto & item : pair.second)
            {
                if( prev == NULL )
                {
                    prev = & item ;
                    continue ;
                }

                (*out)<<prev->base<<'\t'<<item.base<<'\t'
                    <<prev->Rctg<<'\t'
                    <<item.Lctg<<'\t'
                    <<prev->Rquality<<'\t'
                    <<item.Lquality<<'\n';
                prev = & item ;
            }
        }
        delete out ;
    }

    void PrintGapOO()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.gap_oo());
        if( out == NULL )
            FATAL( " failed to open xxx.gap_oo for write ");
        for( const auto & pairs: gaps )
        {
            for( const auto & gap : pairs.second )
            {
                (*out)<<gap.prev<<'\t'<<gap.next<<'\t'
                    <<gap.data.true_prev<<'\t'
                    <<gap.data.true_next<<'\t'
                    <<int (gap.data.value >100 ? 10000 :  (gap.data.value * 100)) 
                    <<'\t'<<int(gap.data.sim * 10000)
                    <<'\n';
            }
        }
        delete out ;
    }
} config;

int main(int argc, char **argv)
{
    //step 0 Parse parmeters...
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, " In xxx.mintree_trunk_linear , xxx.bin_cluster ; xxx.gap_order");
    DEFINE_ARG_OPTIONAL( int , rank , " rank to detect gap ","3");
    END_PARSE_ARGS;

    config.Init( prefix.to_string());
    config.rank = rank.to_int() ;
    config.LoadTrunk();
    config.LoadBinRelationArrayFromFile();
    
    //config.GetGapSim();
    config.BuildSims();
    //config.CalcAll() ;
    config.CalcAll1();
    //config.PrintGapOO() ;
    config.PrintGapOO1() ;
    return 0;
}
