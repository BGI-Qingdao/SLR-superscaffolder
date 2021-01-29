#include "utils/args/argsparser.h"
#include "utils/log/log.h"
#include "utils/log/logfilter.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/string/stringtools.h"
#include "utils/misc/Error.h"
#include "utils/misc/freq.h"
#include "utils/misc/flags.h"

#include "utils/misc/fileName.h"
#include "utils/misc/mapHelper.h"
#include "utils/soap2/contigIndex.h"

#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include "stLFR/ScaffInfo.h"

#include <string.h>

#include "utils/interval/Interval.h"

struct AppConfig
{
    typedef BGIQD::INTERVAL::Interval<float, BGIQD::INTERVAL::IntervalType::Left_Open_Right_Close> SimArea;

    std::map<SimArea, int> gapArea;

    BGIQD::MISC::FileNames fName ;
    BGIQD::LOG::logger loger;

    BGIQD::MISC::Freq<std::string> pfreq;
    BGIQD::MISC::Freq<int>gapFreq;
    int gap_trunk ;
    int gap_petrunk ;
    int gap_pe ;
    bool ptest;
    bool ptest1;

    struct TrueContig
    {
        public:
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
        public:
            float cluster_value ;
            ADD_A_FLAG(4, has_pe);
            ADD_A_FLAG(2, use_pe);
            ADD_A_FLAG(11, 2_pe);
            ADD_A_FLAG(7, pe_confilict);

            ADD_A_FLAG(3, use_sim);
            ADD_A_FLAG(5, has_sim);
            ADD_A_FLAG(13,use_sim_weight);
            ADD_A_FLAG(12,2_sim);
            ADD_A_FLAG(8, sim_confilict);

            ADD_A_FLAG(6, has_pe_sim);
            ADD_A_FLAG(9, pe_sim_confilict);

            ADD_A_FLAG(1, use_basic);
            ADD_A_FLAG(10, no_confilict);
            ADD_A_FLAG(14, use_num);

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
            std::string ToString()  const 
            {
                std::ostringstream ost;
                ost<<basic<<'\t'<<Value()<<'\t'
                    <<sim_prev_ask<<'\t'<<sim_next_ask<<
                    '\t'<<pe_prev_ask<<'\t'<<pe_next_ask;
                return ost.str();
            }
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
                BGIQD::MISC::MapHelper< std::map<unsigned int,int> > Helper;
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
                    f_value = basic ;
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
                        Set_use_pe();
                        return ;
                    }
                    if( Is_2_sim() && !Is_sim_confilict() )
                    {
                        f_value = sim_value() ;
                        Set_use_sim();
                        return ;
                    }
                    auto p1 = vs.begin() ;
                    auto p2 = vs.rbegin() ;
                    if( p1->second != p2->second )
                    {
                        if( p1->second > p2->second )
                            f_value = p1->first ;
                        else
                            f_value = p2->first;
                        Set_use_num();
                        return ;
                    }
                    else
                    {
                        if( Is_has_pe() )
                        {
                            if( ! Is_pe_confilict() )
                            {
                                f_value = pe_value() ;
                                Set_use_pe();
                                return ;
                            }
                        }

                        if( sim_prev_weigth > sim_next_weight )
                        {
                            f_value = sim_prev_ask ;
                        }
                        else
                        {
                            f_value = sim_next_ask ;
                        }
                        Set_use_sim_weight();
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


    BGIQD::SOAP2::ContigIndexMap contigIndexs ;
    void LoadContigIndex()
    {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fName.ContigIndex());
        if(in == NULL)
            FATAL(" failed to open xxx.contigIndex for read!!! ");
        contigIndexs.LoadContigIndexs(*in);
        delete in ;
        contigIndexs.BuildReverseCompleteContigs();
    }
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
    std::vector<ContigOrientation> scaffs;
    int K ;
    void LoadGapArea()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.gap_area());
        if( in == NULL )
            WARN(" failed to open xxx.gap_area for read!!! ");

        if( in == NULL )
        {
            return ;
        }
        std::string line;
        float prev = 1.0f;
        int gap_prev ;
        while( ! std::getline(*in,line).eof() )
        {
            std::istringstream ist(line);
            int gap ; float s ;
            ist>>gap>>s ;
            SimArea tmp(s,prev);
            prev = s ;
            gap_prev = gap ;
            gapArea[tmp] = gap ;
        }
        SimArea tmp(-0.1,prev);
        gapArea[tmp] =  gap_prev;
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load gap_area done "<<BGIQD::LOG::lend() ;
    }

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


    int GetGapLen(float s)
    {
        for( const auto p : gapArea )
        {
            if( p.first.IsContain(s) )
                return p.second ;
        }
        assert(0);
        return -1 ;
    }
    void LoadPEFill()
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.trunk_fill());
        if( in == NULL )
            WARN(" failed to open xxx.trunk_fill for read!!! ");

        if( in == NULL )
            return ;
        auto fill = [this](const std::string & line) ->void
        {
            BGIQD::stLFR::GapFill fill;
            fill.InitFromString(line);
            pefills[fill.prev] = fill ;
            pefills[fill.prev+1] = fill ;
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
                a_scaff[i].cluster_value = float(fill.extra[2])/1000000.0f;
                a_scaff[i].SetValue(fill.true_prev,TrueContig::ValueType::SIM,TrueContig::ValueOrder::NEXT);
                a_scaff[i+1].SetValue(fill.true_next,TrueContig::ValueType::SIM,TrueContig::ValueOrder::PREV);
                a_scaff[i].SetWeight(fill.extra[0],TrueContig::ValueType::SIM,TrueContig::ValueOrder::NEXT);
                a_scaff[i+1].SetWeight(fill.extra[1],TrueContig::ValueType::SIM,TrueContig::ValueOrder::PREV);
                if( pefills.find(gap.prev) != pefills.end() )
                {
                    const auto & pe_fill = pefills[gap.prev];
                    a_scaff[i].SetValue(pe_fill.true_prev,TrueContig::ValueType::PE,TrueContig::ValueOrder::NEXT);
                    a_scaff[i+1].SetValue(pe_fill.true_next,TrueContig::ValueType::PE,TrueContig::ValueOrder::PREV);
                    if( pe_fill.extra.size() > 2 )
                    {
                        a_scaff[i].pe_fill.insert(
                                a_scaff[i].pe_fill.end()
                                , (pe_fill.extra.begin() +1 )
                                , (pe_fill.extra.end() -1 )
                                );
                        pfreq.Touch("With pe step");
                    }
                }
            }

            for( auto & tc : a_scaff )
            {
                tc.CalcValue();
                pfreq.Touch("PointNum");
                if(tc.Is_has_pe() )
                {
                    pfreq.Touch("Has PE");
                }
                if(tc.Is_pe_confilict() )
                {
                    pfreq.Touch("PE confilict");
                }
                if(tc.Is_sim_confilict() )
                {
                    pfreq.Touch("SIM confilict");
                }
                if(tc.Is_pe_sim_confilict() )
                {
                    pfreq.Touch("PE SIM confilict");
                }
                if(tc.Is_has_sim() )
                {
                    pfreq.Touch("Has SIM");
                }
                if(tc.Is_no_confilict() )
                {
                    pfreq.Touch("No confilict");
                }
                if( tc.Is_2_pe() )
                {
                    pfreq.Touch("2 PE");
                }
                if( tc.Is_2_sim() )
                {
                    pfreq.Touch("2 SIM");
                }
                if( tc.Is_use_pe() )
                {
                    pfreq.Touch("USE PE");
                }
                if( tc.Is_use_sim() )
                {
                    pfreq.Touch("USE SIM");
                }
                if( tc.Is_use_sim_weight() )
                {
                    pfreq.Touch("USE SIM WEIGHT");
                }
                if( tc.Is_use_num() )
                {
                    pfreq.Touch("USE NUM");
                }
                if( tc.Is_use_basic() )
                {
                    pfreq.Touch("USE Basic");
                }
            }
            scaffs.push_back(a_scaff);
        }
        loger<<BGIQD::LOG::lstart()<<'\n'<<pfreq.ToString()<<BGIQD::LOG::lend();
    }

   // struct ContigInScaffDetails
   // {
   //     unsigned int base ;
   //     bool plus ;
   //     int gap ;
   // };

    BGIQD::stLFR::ScaffInfoHelper helper ;
    void BuildScaff()
    {
        int id = 0 ;
        for( const auto & a_scaff : scaffs )
        {
            id++ ;
            auto & a_scaffold = helper.all_scaff[id];
            a_scaffold.scaff_id = id ;
            int index = 1 ;
            int pos = 1 ;
            for(size_t i = 0 ; i < a_scaff.size() ; i++ )
            {
                const auto & tc = a_scaff[i];
                unsigned int contig = tc.Value() ;
                if( i == a_scaff.size() - 1 )
                {
                    BGIQD::stLFR::ContigDetail tmp;
                    tmp.contig_id = tc.basic ;
                    tmp.orientation = (tc.basic==contig);
                    tmp.gap_size = 0 ;
                    tmp.contig_len = contigIndexs.GetContigIndex(tc.basic).length ;
                    tmp.scaff_index = index  ;
                    tmp.scaff_id = id ;
                    tmp.start_pos = pos ;
                    a_scaffold.a_scaff.push_back(tmp);
                    break;
                }

                int fill = min_fill ;
                if( tc.pe_fill.empty() )
                {
                    if (tc.pe_next_ask == 0)
                    {
                        if( gapArea.size() > 0 )
                        {
                            fill = GetGapLen(tc.cluster_value) ;
                        }
                        else
                        {
                            fill = gap_trunk ;
                        }
                    }
                    else
                    {
                        fill = gap_petrunk ;
                    }
                }
                else
                {
                    fill = gap_pe ;
                }
                if( fill < min_fill )
                    fill = min_fill ;

                BGIQD::stLFR::ContigDetail tmp;
                tmp.contig_id = tc.basic ;
                tmp.orientation = (tc.basic==contig);
                tmp.gap_size = fill ;
                tmp.contig_len = contigIndexs.GetContigIndex(tc.basic).length ;
                tmp.scaff_index = index  ;
                tmp.scaff_id = id ;
                tmp.start_pos = pos ;
                a_scaffold.a_scaff.push_back(tmp);
                index ++ ;
                pos += (fill+ tmp.contig_len);

                if( ! tc.pe_fill.empty() )
                {
                    for( unsigned int x : tc.pe_fill )
                    {
                        int fill_pe = gap_pe ;
                        if( fill_pe  < min_fill )
                            fill_pe  = min_fill ;
                        unsigned int base = contigIndexs.BaseId(x);
                        BGIQD::stLFR::ContigDetail tmp;
                        tmp.contig_id = base;
                        tmp.orientation = (base == x );
                        tmp.gap_size = fill_pe ;
                        tmp.contig_len = contigIndexs.GetContigIndex(base).length ;
                        tmp.scaff_index = index  ;
                        tmp.scaff_id = id ;
                        tmp.start_pos = pos ;
                        a_scaffold.a_scaff.push_back(tmp);
                        index ++ ;
                        pos += (fill_pe+ tmp.contig_len);
                    }
                }
            }
        }
    }

    std::map< unsigned int , float > to_next_sim ;
    void LoadGapSim()
    {
        auto in  = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fName.gap_sim());
        if( in == NULL )
            WARN(" failed to open xxx.gap_sim for read!!! ");

        if( in == NULL )
        {
            return ;
        }
        std::string line;
        while( ! std::getline(*in,line).eof() )
        {
            int contig , sim ;
            std::istringstream ist(line);
            ist>>contig>>sim ;
            to_next_sim[contig]= float(sim)/1000000.0f ;
        }
    }

    void RewriteSim()
    {
        for( auto & a_scaff : scaffs )
        {
            for( auto & a_contig : a_scaff )
            {
                if( to_next_sim.find( a_contig.basic )
                        ==  to_next_sim.end())
                    continue ;
                a_contig.cluster_value = to_next_sim.at(a_contig.basic) ;
            }
        }
    }

    void PrintScaffInfo()
    {
        auto out1 = BGIQD::FILES::FileWriterFactory::
            GenerateWriterFromFileName(fName.scaff_infos());
        if( out1 == NULL )
            FATAL(" failed to open xxx.scaff_infos to write ");
        helper.PrintAllScaff(*out1);
        delete out1;
    }
    int min_contig ;
    int min_fill ;
} config ;

int main(int argc, char **argv)
{
    //step 0 Parse parmeters...
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , prefix, "prefix of file name\n\
                                                        In\n\
                                                            xxx.mintree_trunk_linear ;\n\
                                                            xxx.contig ;\n\
                                                            xxx.gap_oo ;\n\
                                                            xxx.gap_area;(optional)\n\
                                                            xxx.trunk_fill;(optional)\n\
                                                        Out\n\
                                                            xxx.scaff_infos\n\
                                                        ");
        DEFINE_ARG_OPTIONAL( int , gap_trunk, "gap in trunk" , "5000");
        DEFINE_ARG_OPTIONAL( int , gap_petrunk, "gap in trunk and has pe conn" , "1");
        DEFINE_ARG_OPTIONAL( int , gap_pe, "gap in pe" , "1");
        DEFINE_ARG_OPTIONAL( int , min_gap, "min gap size " , "1");
    END_PARSE_ARGS;

    config.min_fill = min_gap.to_int();
    config.gap_trunk = gap_trunk.to_int();
    config.gap_pe = gap_pe.to_int();
    config.gap_petrunk = gap_petrunk.to_int();
    config.Init(prefix.to_string());

    config.LoadContigIndex();
    config.LoadTrunk();
    config.LoadGapOO();
    config.LoadGapArea();
    config.LoadGapSim();
    config.LoadPEFill();
    config.BuildContigOrientation();
    config.RewriteSim();
    config.BuildScaff();
    config.PrintScaffInfo() ;
    return 0 ;
}
