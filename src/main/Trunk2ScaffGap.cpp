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
#include "soap2/contigIndex.h"

#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <string.h>
#include "algorithm/interval/Interval.h"

#include "biocommon/fasta/fasta.h"

struct AppConfig
{
    typedef BGIQD::INTERVAL::Interval<float, BGIQD::INTERVAL::IntervalType::Left_Open_Right_Close> SimArea;

    typedef BGIQD::FASTA::Fasta<BGIQD::FASTA::ScaffSplitGapHead> GapFasta;

    std::map<SimArea, int> gapArea;

    BGIQD::SOAP2::FileNames fName ;
    BGIQD::LOG::logger loger;

    BGIQD::SOAP2::ContigIndexMap contigIndexs ;
    BGIQD::SOAP2::ContigFastAMap contig_fasta_map;
    BGIQD::FREQ::Freq<std::string> pfreq;
    BGIQD::FREQ::Freq<std::string> gapTypeFreq;
    BGIQD::FREQ::Freq<int>gapFreq;
    int gap_trunk ;
    int gap_petrunk ;
    int gap_pe ;

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

    typedef std::vector<TrueContig> ContigOrientation;

    void LoadContigIndex()
    {
        auto in = BGIQD::FILES::FileReaderFactory::
            GenerateReaderFromFileName(fName.ContigIndex());
        if(in == NULL)
            FATAL(" failed to open xxx.contigIndex for read!!! ");
        contigIndexs.LoadContigIndexs(*in);
        delete in ;
    }

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
        contigMap.K = K;
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

    void BuildScaffGaps()
    {
        auto out = BGIQD::FILES::FileWriterFactory
                   ::GenerateWriterFromFileName(fName.scaff_gap2filler_seqs());
        if( out == NULL )
            FATAL(" failed to open xxx.scaff_gap2filler_seqs to write ");

        int scaff_id = 0 ;

        for( const auto & a_scaff : scaffs )
        {
            scaff_id++ ;
            std::string line ;
            GapFasta tmp;
            int gap_index = 0 ;
            auto add_c1 = [&tmp, &scaff_id ,&gap_index ](
                    unsigned int base ,
                    unsigned int contig ,
                    const BGIQD::SOAP2::ContigFastA &c )
            {
                tmp.Reset();
                tmp.head.scaff_id = scaff_id ;
                gap_index ++ ;
                tmp.head.gap_index = gap_index ;
                tmp.head.prev_contig = contig ;
                tmp.head.prev_base_contig = base;
                tmp.Clean_UnSet();
                tmp.Set_Set_head();
                tmp.AddSeq(c.K);
                tmp.AddSeq(c.linear);
            };

            auto add_ns= [&tmp]( int n , BGIQD::FASTA::ScaffSplitGapHead::GapType  type  )
            {
                tmp.AddSeq(std::string(n,'N'));
                tmp.head.gap_type = type ;
            };

            auto add_c2 = [&tmp, &out](
                    unsigned int base ,
                    unsigned int contig, 
                    const BGIQD::SOAP2::ContigFastA &c )
            {
                tmp.head.next_contig = contig;
                tmp.head.next_base_contig = base;
                tmp.AddSeq(c.K);
                tmp.AddSeq(c.linear);
                (*out)<<tmp.head.Head()<<'\n';
                (*out)<<tmp.seq.Seq(100);
            };

            for( size_t i = 0 ; i < a_scaff.size() ; i++ )
            {
                const auto & tc = a_scaff[i];
                unsigned int contig = tc.Value() ;
                if( i > 0 )
                {
                    add_c2(tc.basic ,contig , contigMap.contigs[contig]);
                }
                if( i != a_scaff.size() - 1 )
                {
                    add_c1(tc.basic, contig,contigMap.contigs[contig]);
                    BGIQD::FASTA::ScaffSplitGapHead::GapType type ;
                    int fill = min_fill ;
                    if( gapArea.size() > 0 )
                    {
                        fill = GetGapLen(tc.cluster_value) ;
                        gapFreq.Touch(fill);
                        // Use default value ;
                        gapTypeFreq.Touch("TrunkGap");
                        type = BGIQD::FASTA::ScaffSplitGapHead::GapType::TRUNK ;
                    }
                    if( tc.pe_fill.empty() )
                    {
                        if (tc.pe_next_ask == 0)
                        {
                            if( gapArea.empty() )
                            {
                                fill = gap_trunk ;
                                type = BGIQD::FASTA::ScaffSplitGapHead::GapType::TRUNK;
                            }
                        }
                        else
                        {
                            gapTypeFreq.Touch("TrunkGap_PE_LINK");
                            fill = gap_petrunk ;
                            type = BGIQD::FASTA::ScaffSplitGapHead::GapType::PE_TRUNK;
                        }
                    }
                    else
                    {
                        fill = gap_pe ;
                        gapTypeFreq.Touch("TrunkGap_PE_FILL");
                        type = BGIQD::FASTA::ScaffSplitGapHead::GapType::PE_TRUNK;
                    }
                    if( fill < min_fill )
                        fill = min_fill ;
                    if( fill <= gap_pe )
                        type = BGIQD::FASTA::ScaffSplitGapHead::GapType::PE_TRUNK;
                    add_ns(fill,type);

                    if( ! tc.pe_fill.empty() )
                    {
                        for( unsigned int x : tc.pe_fill )
                        {
                            add_c2(contigIndexs.BaseId(x),x,contigMap.contigs[x]);
                            add_c1(contigIndexs.BaseId(x),x,contigMap.contigs[x]);
                            int fill_pe = gap_pe ;
                            if( fill_pe < min_fill )
                                fill_pe = min_fill ;
                            add_ns(fill,BGIQD::FASTA::ScaffSplitGapHead::GapType::PE);
                        }
                    }
                }
            }
        }

        delete out;
        loger<<BGIQD::LOG::lstart() << "Build scaff done "<<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart() << " gap freq "
            <<gapFreq.ToString()
            <<BGIQD::LOG::lend() ;

        loger<<BGIQD::LOG::lstart() << " gap type freq "
            <<gapTypeFreq.ToString()
            <<BGIQD::LOG::lend() ;
    }
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
                                                            xxx.scaff_seqs\n\
                                                        ");
        DEFINE_ARG_REQUIRED(int, K, " kvalue ");
        DEFINE_ARG_OPTIONAL( int , gap_trunk, "gap in trunk" , "5000");
        DEFINE_ARG_OPTIONAL( int , gap_petrunk, "gap in trunk and has pe conn" , "300");
        DEFINE_ARG_OPTIONAL( int , gap_pe, "gap in pe" , "50");
        DEFINE_ARG_OPTIONAL( int , min_gap, "min gap size " , "11");
    END_PARSE_ARGS;

    config.min_fill = min_gap.to_int();
    config.Init(prefix.to_string());
    config.K = K.to_int();
    config.gap_trunk = gap_trunk.to_int();
    config.gap_pe = gap_pe.to_int();
    config.gap_petrunk = gap_petrunk.to_int();
    config.LoadContigIndex();
    config.LoadTrunk();
    config.LoadGapOO();
    config.LoadGapArea();
    config.LoadPEFill();
    config.BuildContigOrientation();
    config.LoadContigs();
    config.BuildScaffGaps();

    return 0 ;
}
