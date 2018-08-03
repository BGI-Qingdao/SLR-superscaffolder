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
    END_PARSE_ARGS;

    config.Init( prefix.to_string());
    config.LoadTrunk();
    config.LoadBinRelationArrayFromFile();
    config.GetGapSim();
    config.CalcAll() ;
    config.PrintGapOO() ;
    return 0;
}
