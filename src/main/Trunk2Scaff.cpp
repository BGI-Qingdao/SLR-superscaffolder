#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"
#include "common/stl/mapHelper.h"
#include "common/freq/freq.h"
#include "common/flags/flags.h"

#include "soap2/fileName.h"
#include "soap2/contigFasta.h"
#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <string.h>

struct AppConfig
{
    BGIQD::SOAP2::FileNames fName ;
    BGIQD::LOG::logger loger;

    BGIQD::SOAP2::ContigFastAMap contig_fasta_map;

    struct TrueContig
    {
        private:
            unsigned int sim_prev_ask;
            unsigned int sim_next_ask;

            unsigned int pe_prev_ask;
            unsigned int pe_next_ask;

            unsigned int basic ;
            unsigned int f_value ;

            int sim_prev_weigth ;
            int sim_next_weight ;

            int pe_prev_weigth ;
            int pe_next_weight ;

            FLAGS_INT;
            ADD_A_FLAG(4, has_pe);
            ADD_A_FLAG(2, use_pe);
            ADD_A_FLAG(11, 2_pe);
            ADD_A_FLAG(7, pe_confilict);

            ADD_A_FLAG(3, use_sim);
            ADD_A_FLAG(5, has_sim);
            ADD_A_FLAG(12,2_sim);
            ADD_A_FLAG(8, sim_confilict);

            ADD_A_FLAG(6, has_pe_sim);
            ADD_A_FLAG(9, pe_sim_confilict);

            ADD_A_FLAG(1, use_basic);
            ADD_A_FLAG(10, no_confilict);

            unsigned int pe_value() const 
            {
                if( pe_prev_ask != 0 )
                    return pe_prev_ask ;
                if( pe_next_ask != 0 )
                    return pe_next_ask ;
                assert(0);
                return 0;
            }
            unsigned int sim_value() const 
            {
                if( sim_prev_ask != 0 )
                    return sim_prev_ask ;
                if( sim_next_ask != 0 )
                    return sim_next_ask ;
                assert(0);
                return 0;
            }
        public :
            std::vector<unsigned int> pe_fill ;

            void Init(unsigned int bs)
            {
                basic = bs ;
                pe_next_ask = 0 ;
                pe_prev_ask = 0 ;
                sim_next_ask = 0 ;
                sim_prev_ask = 0 ;
                pe_prev_weigth = 0 ;
                pe_next_weight = 0 ;
                sim_prev_weigth = 0 ;
                sim_next_weight = 0 ;
                f_value = 0 ;
            }

            enum ValueType
            {
                Basic,
                PE,
                SIM,
            };
            enum ValueOrder
            {
                PREV,
                NEXT
            };

            void SetWeight( unsigned int value , ValueType t , ValueOrder o )
            {
                if( t == ValueType::PE && o == ValueOrder::PREV )
                    pe_prev_weigth = value ;
                if( t == ValueType::PE && o == ValueOrder::NEXT )
                    pe_next_weight = value ;

                if( t == ValueType::SIM && o == ValueOrder::PREV )
                    sim_prev_weigth = value ;
                if( t == ValueType::SIM && o == ValueOrder::NEXT )
                    sim_next_weight = value ;
            }

            void SetValue(unsigned int value , ValueType t , ValueOrder o )
            {
                if( t == ValueType::PE && o == ValueOrder::PREV )
                    pe_prev_ask = value ;
                if( t == ValueType::PE && o == ValueOrder::NEXT )
                    pe_next_ask = value ;

                if( t == ValueType::SIM && o == ValueOrder::PREV )
                    sim_prev_ask = value ;
                if( t == ValueType::SIM && o == ValueOrder::NEXT )
                    sim_next_ask = value ;
            }

            bool HasValue() const { return pe_prev_ask >0 
                || pe_next_ask > 0
                    || sim_next_ask > 0
                    || sim_prev_ask > 0
                    ; }

            void CalcValue()
            {
                std::map<unsigned int,int> vs ;
                int pe = 0 ;
                int sim = 0 ;
                BGIQD::STL::MapHelper< std::map<unsigned int,int> > Helper;
                if( pe_next_ask != 0 )
                {
                    Set_has_pe();
                    Helper.Incr(vs,pe_next_ask,1);
                    pe++ ;
                }
                if( pe_prev_ask != 0 )
                {
                    Set_has_pe();
                    Helper.Incr(vs,pe_prev_ask,1);
                    pe++ ;
                }
                if( sim_next_ask != 0 )
                {
                    Set_has_sim();
                    Helper.Incr(vs,sim_next_ask,1);
                    sim++ ;
                }
                if( sim_prev_ask != 0 )
                {
                    Set_has_sim() ;
                    Helper.Incr(vs,sim_prev_ask,1);
                    sim ++ ;
                }
                if( pe == 2 )
                {
                    if( pe_next_ask != pe_prev_ask )
                        Set_pe_confilict() ;
                    Set_2_pe();
                }
                if( sim == 2 )
                {
                    if (sim_next_ask != sim_prev_ask )
                        Set_sim_confilict() ;
                    Set_2_sim();
                }
                if( ( Is_has_pe() && ! Is_pe_confilict() ) && ( Is_has_sim() && ! Is_sim_confilict() ) )
                {
                    if( pe_value() != sim_value() )
                    {
                        Set_pe_sim_confilict() ;
                    }
                }
                if (vs.size() == 0 )
                {
                    Set_use_basic();
                    f_value = 0 ;
                    return ;
                }
                else if( vs.size() == 1 )
                {
                    Set_no_confilict();
                    f_value = vs.begin()->first ;
                    return ;
                }
                else
                {
                    if( Is_2_pe() && !Is_pe_confilict() )
                    {
                        f_value = pe_value() ;
                        return ;
                    }
                    if( Is_2_sim() && !Is_sim_confilict() )
                    {
                        f_value = sim_value() ;
                        return ;
                    }
                    auto p1 = vs.begin() ;
                    auto p2 = vs.rbegin() ;
                    if( p1->second != p2->second )
                    {
                        if( p1->second > p2->second )
                            f_value = p1->first ;
                        else
                            f_value = p2->second ;
                        return ;
                    }
                    else
                    {
                        if( Is_has_pe() )
                        {
                            if( Is_pe_confilict() )
                            {
                                f_value = pe_value() ;
                                return ;
                            }
                        }
                        else
                        { // PE 1 + 1^ , SIM 1 + 1^
                            if( sim_prev_weigth > sim_next_weight )
                            {
                                f_value = sim_prev_ask ;
                            }
                            else
                            {
                                f_value = sim_next_ask ;
                            }
                        }
                    }
                }
                return ;
            }
            unsigned int Value() const 
            {
                assert( f_value > 0 && (f_value == basic || f_value == basic+1));
                return f_value;
            }
    };

    typedef std::vector<TrueContig> ContigOrientation;
    // struct GapPos
    // {
    //     int trunkId ;
    //     int gapId;
    // };
    struct GapExtra
    {
        unsigned int true_prev ;
        unsigned int true_next ; 
        bool is_PE ;
        int value ;
    };

    std::map<int , std::vector<BGIQD::stLFR::TrunkGap<GapExtra> >> gaps;
    //  std::map<unsigned int , GapPos> contigPos;
    std::map<unsigned int, BGIQD::stLFR::GapFill> gapfills;
    std::map<unsigned int, BGIQD::stLFR::GapFill> pefills;
    BGIQD::SOAP2::ContigFastAMap contigMap;
    std::vector<ContigOrientation> scaffs;

    void Init(const std::string & prefix)
    {
        fName.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("Trunk2Scaff ",BGIQD::LOG::loglevel::INFO, loger);
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

    //    void BuildGapPos()
    //    {
    //        for( const auto & p : gaps )
    //        {
    //            for( size_t i = 0 ; i < p.second.size() ; i++)
    //            {
    //                const auto gap = p.second[i];
    //                contigPos[gap.prev] = GapPos{ p.first , (int)i } ;
    //            }
    //        }
    //    }

    void LoadGapOO()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.gap_oo());
        if( in == NULL )
            FATAL(" failed to open xxx.gap_oo for read!!! ");
        auto fill = [this](const std::string & line) ->void
        {
            BGIQD::stLFR::GapFill fill;
            fill.InitFromString(line);
            gapfills[fill.prev] = fill ;
            gapfills[fill.prev+1] = fill ;
        };

        BGIQD::FILES::FileReaderFactory::EachLine(*in,fill);
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load GapOO done "<<BGIQD::LOG::lend() ;
    }

    void LoadContigs()
    {
        contigMap.LoadContig(fName.contig());
        contigMap.buildCompeleReverse();
        loger<<BGIQD::LOG::lstart() << "Load Contig done "<<BGIQD::LOG::lend() ;
    }

    /*
    void BuildScaff()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.scaff_seqs());
        if( out == NULL )
            FATAL(" failed to open xxx.scaff_seqs to write ");

        std::vector<std::vector<unsigned int> >scaffs ;

        for( const auto & pair : gaps )
        {
            std::vector<unsigned int> a_scaff;
            for(size_t i = 0 ; i< pair.second.size() ; i++ )
            {
                const auto & gap = pair.second[i];
                const auto & fill = gapfills[gap.prev];

                if( a_scaff.empty() )
                {
                    a_scaff.push_back(fill.true_prev); 
                    a_scaff.push_back(fill.true_next); 
                }
                else
                {
                    if( *(a_scaff.rbegin()) == fill.true_prev ) 
                    {
                        a_scaff.push_back(fill.true_next);
                    }
                    else
                    {
                        auto & fill_prev = gapfills[*(a_scaff.rbegin()+1)] ;
                        int value = *fill_prev.extra.begin();
                        int value1 =* fill.extra.begin() ;
                        if( value1 > value )
                        {
                            *(a_scaff.rbegin()) = fill.true_prev ;
                        }
                        a_scaff.push_back(fill.true_next);
                    }
                }
            }
            scaffs.push_back(a_scaff);
        }
        int id = 0 ;
        for( const auto & a_scaff : scaffs )
        {
            id++ ;
            std::string line ;
            (*out)<<">scaffold"<<id<<"\t20.5\n";
            for(size_t i = 0 ; i < a_scaff.size() ; i++ )
            {
                unsigned int contig = a_scaff[i] ;
                line+=contigMap.contigs[contig].K;
                line+=contigMap.contigs[contig].linear;
                if( i != a_scaff.size() - 1 )
                {
                    line += std::string(5000,'N');
                }
            }
            for( int i = 0 ; i < (int)line.size() ; i++ )
            {
                (*out)<<line[i];
                if( i % 100 == 0 || i ==(int) line.size() -1 )
                {
                    (*out)<<'\n';
                }
            }
        }
        delete out;
        loger<<BGIQD::LOG::lstart() << "Build scaff done "<<BGIQD::LOG::lend() ;
    }
    */
    void LoadPEFill()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.trunk_fill());
        if( in == NULL )
            FATAL(" failed to open xxx.trunk_fill for read!!! ");
        auto fill = [this](const std::string & line) ->void
        {
            BGIQD::stLFR::GapFill fill;
            fill.InitFromString(line);
            gapfills[fill.prev] = fill ;
            gapfills[fill.prev+1] = fill ;
        };

        BGIQD::FILES::FileReaderFactory::EachLine(*in,fill);
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load PE Trunk fill done "<<BGIQD::LOG::lend() ;
    }

    void BuildContigOrientation()
    {
        for( const auto & pair : gaps )
        {
            ContigOrientation a_scaff ;
            for(size_t i = 0 ; i< pair.second.size() ; i++ )
            {
                TrueContig tc;
                const auto & gap = pair.second[i];
                if( i == 0 )
                {
                    TrueContig tc0;
                    tc0.Init(gap.prev);
                    a_scaff.push_back(tc0);
                }
                tc.Init(gap.next);
                a_scaff.push_back(tc);
            }
            for(size_t i = 0 ; i< pair.second.size() ; i++ )
            {
                const auto & gap = pair.second[i];
                const auto & fill = gapfills[gap.prev];
                a_scaff[i].SetValue(fill.true_prev,TrueContig::ValueType::SIM,TrueContig::ValueOrder::NEXT);
                a_scaff[i+1].SetValue(fill.true_next,TrueContig::ValueType::SIM,TrueContig::ValueOrder::PREV);
                a_scaff[i].SetWeight(*fill.extra.begin(),TrueContig::ValueType::SIM,TrueContig::ValueOrder::NEXT);
                a_scaff[i+1].SetWeight(*fill.extra.begin(),TrueContig::ValueType::SIM,TrueContig::ValueOrder::PREV);
                if( pefills.find(gap.prev) != pefills.end() )
                {
                    const auto & pe_fill = pefills[gap.prev];
                    a_scaff[i].SetValue(pe_fill.true_prev,TrueContig::ValueType::SIM,TrueContig::ValueOrder::NEXT);
                    a_scaff[i+1].SetValue(pe_fill.true_next,TrueContig::ValueType::SIM,TrueContig::ValueOrder::PREV);
                    if( pe_fill.extra.size() > 2 )
                    {
                        a_scaff[i].pe_fill.insert(
                                a_scaff[i].pe_fill.end()
                                , (pe_fill.extra.begin() +1 )
                                , (pe_fill.extra.end() -1 )
                                );
                    }
                }
            }
            for( auto & tc : a_scaff )
            {
                tc.CalcValue();
            }
            scaffs.push_back(a_scaff);
        }
    }

    void BuildScaff()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.scaff_seqs());
        if( out == NULL )
            FATAL(" failed to open xxx.scaff_seqs to write ");

        int id = 0 ;
        for( const auto & a_scaff : scaffs )
        {
            id++ ;
            std::string line ;
            (*out)<<">scaffold"<<id<<"\t20.5\n";
            for(size_t i = 0 ; i < a_scaff.size() ; i++ )
            {
                const auto & tc = a_scaff[i];
                unsigned int contig = tc.Value() ;

                line+=contigMap.contigs[contig].K;
                line+=contigMap.contigs[contig].linear;
                if( tc.pe_fill.empty() )
                {
                    if( i != a_scaff.size() - 1 )
                    {
                        line += std::string(5000,'N');
                    }
                }
                else
                {
                    line += std::string(10,'N');
                    for( unsigned int x : tc.pe_fill )
                    {
                        line+=contigMap.contigs[x].K;
                        line+=contigMap.contigs[x].linear;
                        line += std::string(10,'N');
                    }
                }
            }
            for( int i = 0 ; i < (int)line.size() ; i++ )
            {
                (*out)<<line[i];
                if( i % 100 == 0 || i ==(int) line.size() -1 )
                {
                    (*out)<<'\n';
                }
            }
        }
        delete out;
        loger<<BGIQD::LOG::lstart() << "Build scaff done "<<BGIQD::LOG::lend() ;
    }
} config ;

int main(int argc, char **argv)
{
    //step 0 Parse parmeters...
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix, " In xxx.mintree_trunk_linear , xxx.bin_cluster ; xxx.gap_order");
    END_PARSE_ARGS;

    config.Init(prefix.to_string());
    config.LoadTrunk();
    config.LoadGapOO();
    config.LoadPEFill();
    config.BuildContigOrientation();
    config.LoadContigs();
    config.BuildScaff();

    return 0 ;
}
