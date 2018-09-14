#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/flags/flags.h"
#include "common/freq/freq.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "algorithm/multi_key_hash/MultiKeyHash.h"

#include "stLFR/EasySam.h"
#include "stLFR/CBB.h"

#include "soap2/contigGraph.h"
#include "soap2/fileName.h"

#include <algorithm>

struct AppConfig {

    BGIQD::LOG::logger loger;

    int K;

    int min_cout ;

    float min_factor;

    int bmin ;

    float bfac ;

    typedef BGIQD::MultiKeyMap::BiKeyHash<unsigned int , int> PECache;

    PECache pe_cache ;

    std::map<unsigned int , BGIQD::stLFR::ContigBarcodeInfo> cbs;

    void TouchPECahce(unsigned int c1 , unsigned int c2)
    {
        if( pe_cache.Contain(c1,c2) )
        {
            pe_cache.Set(c1 , c2 , pe_cache.At(c1,c2)+1);
        }
        else
        {
            pe_cache.Set(c1 , c2 , 1);
        }
    }

    void TouchPECahce(unsigned int c1 , unsigned int c2, int i)
    {
        pe_cache.Set(c1 , c2 , i);
    }

    BGIQD::SOAP2::GraphEA graph_ea ;
    BGIQD::SOAP2::FileNames fNames;

    struct SeedExternInfo
    {
        BGIQD::stLFR::ContigIndex c_info;
        std::vector<unsigned int> left_extern_path;
        std::vector<unsigned int> right_extern_path;

        std::set<unsigned int> l_group ;
        std::set<unsigned int> r_group ;


        std::set<int> barcodes ;


        int ShareBarcodeCount(const BGIQD::stLFR::ContigBarcodeInfo & info)
        {
            int ret = 0 ;
            for( const auto & pair : info.barcodesOnPos )
            {
                for( int i : pair.second ) 
                    if( barcodes.find(i) != barcodes.end() )
                        ret ++ ;
            }
            return ret ;
        }

        void UpdateBarcodes(const BGIQD::stLFR::ContigBarcodeInfo & info)
        {
            for( const auto & pair : info.barcodesOnPos )
            {
                for( int i : pair.second ) 
                    if( barcodes.find(i) == barcodes.end() )
                        barcodes.insert(i);
            }
        }

        bool AddRightPath( unsigned int c , unsigned int c1 )
        {
            if( r_group.find( c ) != r_group.end()
                    || r_group.find(c1) != r_group.end()) 
                return false ;
            right_extern_path.push_back(c);
            r_group.insert(c);
            l_group.insert(c1);
            return true ;
        }

        bool AddLeftPath( unsigned int c , unsigned int c1)
        {
            if( l_group.find( c ) != l_group.end()
                    || l_group.find(c1) != l_group.end()) 
                return false ;
            left_extern_path.push_back(c);
            l_group.insert(c);
            r_group.insert(c1);
            return true ;
        }

        bool operator < ( const SeedExternInfo & a ) const 
        {
            return c_info.length < a.c_info.length ;
        }

        FLAGS_INT
        ADD_A_FLAG(2,LEndByNoInfo);
        ADD_A_FLAG(3,LEndByConfuse);
        ADD_A_FLAG(4,LEndByNoMoreArc);
        ADD_A_FLAG(5,LEndBy1ArcFork);
        ADD_A_FLAG(6,LEndByNo1Arc);
        ADD_A_FLAG(7,LEndByCircle);
        ADD_A_FLAG(8,LEndByNoBarcode);
        ADD_A_FLAG(9,LEndByBarcodeFailed);

        ADD_A_FLAG(12,REndByNoInfo);
        ADD_A_FLAG(13,REndByConfuse);
        ADD_A_FLAG(14,REndByNoMoreArc);
        ADD_A_FLAG(15,REndBy1ArcFork);
        ADD_A_FLAG(16,REndByNo1Arc);
        ADD_A_FLAG(17,REndByCircle);
        ADD_A_FLAG(18,REndByNoBarcode);
        ADD_A_FLAG(19,REndByBarcodeFailed);

        ADD_A_FLAG(1,UsedAsFill);
        ADD_A_FLAG(21,UsedAsCenter);
        ADD_A_FLAG(22, Has_Extern);
        ADD_A_FLAG(23, Has_Extern_Unless);

        void Init()
        {
            flags = 0 ;
            r_group.insert(c_info.contig);
            l_group.insert(c_info.contig+1);
        }
    };

    std::vector<SeedExternInfo> seeds;
    std::map<unsigned int, unsigned int> seed_ids ;
    std::map<unsigned int, int> seed_index;
    //PE

    void LoadPECache()
    {
        BGIQD::LOG::timer t(loger,"LoadPECache");
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.pe_pairs()) ;
        if( in == NULL )
            FATAL( "open .pe_pairs file to read failed !!! " );

        auto eachline = [this](const std::string & line) ->void 
        {
            BGIQD::EASY_SAM::PEInfo tmp ;
            tmp.InitFromString(line);

            const auto & c1 = graph_ea.edge_array[tmp.contig1];
            const auto & c2 = graph_ea.edge_array[tmp.contig2];
            if( tmp.match_reverse1 ^ tmp.match_reverse2 )
            {
                TouchPECahce(c1.id , c2.id ) ;  
                TouchPECahce(c1.bal_id , c2.bal_id ) ;  
            }
            else
            {
                TouchPECahce(c1.id , c2.bal_id ) ;  
                TouchPECahce(c1.bal_id , c2.id ) ;  
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,eachline);
        delete in ;
    }

    void LoadPECache_1()
    {
        BGIQD::LOG::timer t(loger,"LoadPECache");
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(fNames.contig_pe_conns()) ;
        if( in == NULL )
            FATAL( "open .contig_pe_conns file to read failed !!! " );

        auto eachline = [this](const std::string & line) ->void 
        {
            unsigned int cc1 ,cc2;
            int count ;
            char pos ;
            std::istringstream ist(line);
            ist>>cc1>>cc2>>count>>pos;
            const auto & c1 = graph_ea.edge_array[cc1];
            const auto & c2 = graph_ea.edge_array[cc2];
            if( pos == '+')
            {
                TouchPECahce(cc1,cc2,count);
                TouchPECahce(c1.bal_id,c2.bal_id,count);
            }
            else
            {
                TouchPECahce(c1.bal_id,cc2,count);
                TouchPECahce(cc1,c2.bal_id,count);
            }
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,eachline);
        delete in ;
    }

    void LoadGraphEA()
    {
        BGIQD::LOG::timer t(loger,"LoadGraphEA");
        graph_ea.LoadEdge(fNames.updatedEdge(),K);
        graph_ea.LoadArc(fNames.Arc());
    }

    void LoadSeeds()
    {
        BGIQD::LOG::timer t(loger,"LoadSeeds");
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fNames.pe_seeds()) ;
        if( in == NULL )
            FATAL( "open .pe_seeds file to read failed !!! " );

        std::string line ;
        while( in && !std::getline(*in, line).eof() )
        {
            SeedExternInfo tmp ;
            tmp.c_info.InitFromString(line);
            tmp.Init() ;
            seeds.push_back(tmp);
            seed_ids[tmp.c_info.contig] = tmp.c_info.contig;
            seed_ids[tmp.c_info.contig+1] = tmp.c_info.contig;
        }
        delete in ;
        std::sort(seeds.rbegin() , seeds.rend());
        for( int i = 0; i<(int)seeds.size() ; i++)
        {
            seed_index[seeds[i].c_info.contig] = i ;
        }
    }

    void ExternSeeds(SeedExternInfo & seed)
    {
        struct ChooseData
        {
            int c , c1 ; 
            BGIQD::SOAP2::Arc * a, *a1;
        };

        auto GetBarcodeInfo = [this](unsigned int s )
            -> BGIQD::stLFR::ContigBarcodeInfo &
            {
                const auto & edge = graph_ea.edge_array[s] ;
                if( edge.IsBase() )
                    return cbs[s] ;
                else
                    return cbs[s-1] ;
            };

        auto extractPEData = [&](BGIQD::SOAP2::Arc * arc ) 
        {
            std::vector<std::tuple<int, BGIQD::SOAP2::Arc *> > data;
            do{
                unsigned int to = arc->to ;
                int count = 0 ;
                for( const auto i : seed.r_group )
                {
                    if( pe_cache.Contain( i , to ) )
                    {
                        count += pe_cache.At( i , to );
                    }
                }
                data.push_back( std::make_tuple( count , arc ) );
                arc = arc->next ;
            }while(arc != NULL );
            assert(data.size() >1);
            std::sort(data.rbegin() , data.rend());
            ChooseData ret  ;
            std::tie(ret.c,ret.a) = data[0];
            std::tie(ret.c1,ret.a1) = data[1];
            return ret;
        };

        auto extractBarcodeData = [&]( BGIQD::SOAP2::Arc * arc ) 
        {
            std::vector<std::tuple<int, BGIQD::SOAP2::Arc *> > data;
            do{
                unsigned int to = arc->to ;
                int count =  seed.ShareBarcodeCount(GetBarcodeInfo(to));
                data.push_back( std::make_tuple( count , arc ) );
                arc = arc->next ;
            }while(arc != NULL );
            assert(data.size() >1);
            std::sort(data.rbegin() , data.rend());
            ChooseData ret  ;
            std::tie(ret.c,ret.a) = data[0];
            std::tie(ret.c1,ret.a1) = data[1];
            return ret;
        };

        auto CheckSeeds = [this] ( unsigned int s ) -> void
        {
            if( seed_ids.find(s) != seed_ids.end() )
            {
                unsigned int seed_id = seed_ids.at(s) ;
                int index = seed_index.at(seed_id);
                auto & seed  = seeds[index];
                if( ! seed.Is_UsedAsFill() )
                {
                    seed.Set_UsedAsFill() ;
                }
                else
                {
                    //TODO
                }
            }
        };

        auto & edge = graph_ea.edge_array[seed.c_info.contig];
        seed.Set_UsedAsCenter();
        seed.UpdateBarcodes(cbs[seed.c_info.contig]);
        // Downstream
        do{
            BGIQD::SOAP2::Arc * next_cross = edge.arc ;
            if( next_cross == NULL )
            {
                seed.Set_REndByNo1Arc();
                break;
            }
            while( next_cross != NULL )
            {
                if( next_cross->next == NULL )
                {
                    unsigned int next = next_cross->to ;
                    auto & next_edge = graph_ea.edge_array[next] ;
                    if(! seed.AddRightPath(next_edge.id , next_edge.bal_id) )
                    {
                        seed.Set_REndByCircle() ;
                        break ;
                    }
                    seed.UpdateBarcodes(GetBarcodeInfo(next));
                    CheckSeeds(next);
                    next_cross = next_edge.arc ;
                }
                else
                {
                    BGIQD::SOAP2::Arc * arc = next_cross ;
                    auto ret = extractPEData(arc); 
                    if( ret.c < min_cout -1 || float(ret.c)/float(ret.c1) < min_factor )
                    {
                        auto ret1 = extractBarcodeData(arc);
                        if( ret1.c < bmin || float(ret1.c)/float(ret1.c1) < bfac )
                        {
                            seed.Set_REndByBarcodeFailed();
                            break ;
                        }
                        ret.a = ret1.a ;
                    }
                    auto & next_edge = graph_ea.edge_array[ret.a->to] ;
                    CheckSeeds(ret.a->to);
                    seed.UpdateBarcodes(GetBarcodeInfo(ret.a->to));
                    if(! seed.AddRightPath(next_edge.id , next_edge.bal_id) )
                    {
                        seed.Set_REndByCircle() ;
                        break ;
                    }
                    next_cross = next_edge.arc ;
                }
            }
            if( next_cross == NULL )
            {
                seed.Set_REndByNoMoreArc();
            }
        }while(0);

        // Upstream
        auto & bal_edge = graph_ea.edge_array[edge.bal_id];
        do{
            BGIQD::SOAP2::Arc * next_cross = bal_edge.arc ;
            if( next_cross == NULL )
            {
                seed.Set_LEndByNo1Arc();
                break;
            }
            while( next_cross != NULL )
            {
                if( next_cross->next == NULL )
                {
                    unsigned int next = next_cross->to ;
                    auto & next_edge = graph_ea.edge_array[next] ;
                    if( ! seed.AddLeftPath(next_edge.id , next_edge.bal_id) )
                    {
                        seed.Set_LEndByCircle() ;
                        break ;
                    }
                    CheckSeeds(next);
                    seed.UpdateBarcodes(GetBarcodeInfo(next));
                    next_cross = next_edge.arc ;
                }
                else
                {
                    BGIQD::SOAP2::Arc * arc = next_cross ;
                    auto ret = extractPEData(arc); 
                    if( ret.c < min_cout -1 || float(ret.c)/float(ret.c1) < min_factor )
                    {
                        auto ret1 = extractBarcodeData(arc);
                        if( ret1.c < bmin || float(ret1.c)/float(ret1.c1) < bfac )
                        {
                            seed.Set_LEndByBarcodeFailed();
                            break ;
                        }
                        ret.a = ret1.a ;
                    }
                    auto & next_edge = graph_ea.edge_array[ret.a->to] ;
                    if( ! seed.AddLeftPath(next_edge.id , next_edge.bal_id)) 
                    {
                        seed.Set_LEndByCircle() ;
                        break ;
                    }
                    CheckSeeds(ret.a->to);
                    seed.UpdateBarcodes(GetBarcodeInfo(ret.a->to));
                    next_cross = next_edge.arc ;
                }
            }
            if( next_cross == NULL )
            {
                seed.Set_LEndByNoMoreArc();
            }
        }while(0);
        if( seed.left_extern_path.size() > 0 || seed.right_extern_path.size() > 0 )
            seed.Set_Has_Extern();
    }

    void ExternAll()
    {
        int count = 0 ;
        for( auto & a_seed : seeds )
        {
            count ++ ;
            if( count > 1000 && count % 1000 == 1 )
            {
                loger<<BGIQD::LOG::lstart()
                    <<"Process "<<count<<" seeds now ... "
                    <<BGIQD::LOG::lend();
            }
            if( a_seed.Is_UsedAsFill() )
                continue ;
            ExternSeeds(a_seed);
        }
    }

    void PrintResult()
    {
        auto seed2line = [this]( const SeedExternInfo & seed)
        {
            std::vector<unsigned int> line;
            for( unsigned int i : seed.left_extern_path )
            {
                const auto & edge = graph_ea.edge_array[i];
                line.push_back(edge.bal_id);
            }
            std::vector<unsigned int> ret;
            ret.insert(ret.end() , line.rbegin(), line.rend());
            ret.push_back(seed.c_info.contig);
            ret.insert(ret.end() , seed.right_extern_path.begin() , seed.right_extern_path.end());
            return ret;
        };
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fNames.seed_extern_fill());
        if( out == NULL )
            FATAL( "failed to open xxx.seed_extern_fill to write ");
        for( auto & a_seed : seeds )
        {
            if(!  a_seed.Is_UsedAsCenter() || ! a_seed.Is_Has_Extern())
                continue ;
            if( a_seed.left_extern_path.size() < 2 
                    && a_seed.right_extern_path.size() < 2 )
            {
                a_seed.Set_Has_Extern_Unless() ;
                continue ;
            }
            auto ret = seed2line(a_seed);
            *out<<"l";//ret.size();
            for( unsigned int i : ret )
            {
                *out<<'\t'<<i;
            }
            *out<<'\n';
        }
        delete out ;
    }

    void PrintStatistics()
    {
        typedef BGIQD::FREQ::Freq<std::string> TagFreq;
        TagFreq seedType ;
        TagFreq seedExtType;
        TagFreq seedLEnd ;
        TagFreq seedREnd ;
        for( auto & a_seed : seeds )
        {
            if(! a_seed.Is_UsedAsCenter() ) 
            {
                assert( a_seed.Is_UsedAsFill());
                seedType.Touch("UsedAsFill");
            }
            else
            {
                seedType.Touch("UsedAsCenter");
                if( ! a_seed.Is_Has_Extern() )
                {
                    seedExtType.Touch("FillNone");
                }
                else
                {
                    if( a_seed.Is_Has_Extern_Unless() )
                        seedExtType.Touch("Unless");
                    else
                        seedExtType.Touch("Fill");
                }

                if( a_seed.Is_LEndByConfuse() )
                    seedLEnd.Touch("LEndByConfuse");
                else if ( a_seed.Is_LEndByNo1Arc() )
                    seedLEnd.Touch("LEndByNo1Arc");
                else if ( a_seed.Is_LEndByNoMoreArc() )
                    seedLEnd.Touch("LEndByNoMoreArc");
                else if ( a_seed.Is_LEndByNoInfo() )
                    seedLEnd.Touch("LEndByNoInfo");
                else if ( a_seed.Is_LEndBy1ArcFork() )
                    seedLEnd.Touch("LEndBy1ArcFork");
                else if ( a_seed.Is_LEndByCircle() )
                    seedLEnd.Touch("LEndByCircle");
                else
                    assert(0);
                if( a_seed.Is_REndByConfuse() )
                    seedREnd.Touch("REndByConfuse");
                else if ( a_seed.Is_REndByNo1Arc() )
                    seedREnd.Touch("REndByNo1Arc");
                else if ( a_seed.Is_REndByNoMoreArc() )
                    seedREnd.Touch("REndByNoMoreArc");
                else if ( a_seed.Is_REndByNoInfo() )
                    seedREnd.Touch("REndByNoInfo");
                else if ( a_seed.Is_REndBy1ArcFork() )
                    seedREnd.Touch("REndBy1ArcFork");
                else if ( a_seed.Is_REndByCircle() )
                    seedREnd.Touch("REndByCircle");
                else 
                    assert(0);
            }
        }
        loger<<BGIQD::LOG::lstart()
            <<"seed Used As  Type\n"
            <<seedType.ToString()
            <<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart()
            <<"seed Extern As  Type\n"
            <<seedExtType.ToString()
            <<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart()
            <<"seed left end As  Type\n"
            <<seedLEnd.ToString()
            <<BGIQD::LOG::lend() ;
        loger<<BGIQD::LOG::lstart()
            <<"seed right end As  Type\n"
            <<seedREnd.ToString()
            <<BGIQD::LOG::lend() ;
    }

    void Init( const std::string & prefix )
    {
        fNames.Init(prefix);
        BGIQD::LOG::logfilter::singleton().get("PEGraph",BGIQD::LOG::loglevel::INFO, loger);
        BGIQD::stLFR::ContigIndex::K = K ;
    }

    void LoadBarcodeAtContig()
    {
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName
            (fNames.barcode_at_contig_v1());
        if( in  == NULL )
            FATAL( "xxx.barcode_at_contig_v1 open for read failed !!! ");

        auto parseline = [&](const std::string & line )
        {
            BGIQD::stLFR::ContigBarcodeInfo tmp ;
            tmp.InitFromString(line);
            cbs[tmp.contig_id] = tmp;
        };

        BGIQD::FILES::FileReaderFactory::EachLine(*in , parseline);

        delete in ;
    }

} config ;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix , "prefix of file name \n\
                                                    Input \n\
                                                            xxx.pe_seeds\n\
                                                            xxx.seed_extern_fill\n\
                                                            xxx.updated.edge\n\
                                                            xxx.Arc\n\
                                                    Output \n\
                                                            xxx.seed_extern_fill");
    DEFINE_ARG_REQUIRED(int , kvalue , "kvalue used by SOAP");
    DEFINE_ARG_OPTIONAL(int , min_count, "min PE conn count", "20");
    DEFINE_ARG_OPTIONAL(float ,min_factor, "min factor that biggest conn count / second biggest","3");
    DEFINE_ARG_OPTIONAL(int , min_bcount, "min barcode conn count", "2");
    DEFINE_ARG_OPTIONAL(float ,min_bfactor, "min factor that biggest barcode count / second biggest","2");
    END_PARSE_ARGS

    config.K = kvalue.to_int();
    config.min_cout = min_count.to_int();
    config.min_factor = min_factor.to_float();
    config.bmin = min_bcount.to_int() ;
    config.bfac = min_bfactor.to_float() ;
    config.Init(prefix.to_string());

    config.LoadSeeds();
    config.LoadGraphEA();
    config.LoadPECache_1();
    config.LoadBarcodeAtContig();
    config.ExternAll();
    config.PrintResult();
    config.PrintStatistics();
}
