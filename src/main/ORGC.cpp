#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"

#include "stLFR/ScaffInfo.h"
#include "stLFR/SmithWaterman_OO.h"

#include <string>
#include <vector>
#include <sstream>
#include <cassert>
#include <map>
#include <set>
#include <algorithm>


struct ReadInfo
{
    std::string order;
    char orientation ;
    int pos ;
    std::string ref_id ;
    //bool operator < ( const ReadInfo & o )
    //{
    //    if( read_id != o.read_id )
    //        return read_id< o.read_id ;
    //    else 
    //        return pos < o.pos ;
    //}
};

typedef BGIQD::stLFR::SmithWaterman_OO<ReadInfo> TheSmithWaterman;

struct MapperInfo
{
    std::string ref_id ;
    std::vector<ReadInfo>  read_infos ;

    void InitFromStr( const std::string & line )
    {
        std::istringstream ist(line);
        ist>>ref_id ;
        ReadInfo tmp ;
        while( ! ist.eof() )
        {
            ist>>tmp.order>>tmp.orientation>>tmp.pos;
            read_infos.push_back(tmp);
            reads.insert(tmp.order);
        }
    }

    static std::vector<ReadInfo> CompleteReverse( const std::vector<ReadInfo> & c1 )
    {
        auto change_oo =  [](ReadInfo & i)
        { 
            if(i.orientation == '+') 
                i.orientation = '-' ; 
            else 
                i.orientation = '+' ; 
        };
        std::vector<ReadInfo> ret ; 
        ret.resize(c1.size());
        std::copy_backward(c1.begin(),c1.end(),ret.begin());
        std::for_each(ret.begin() ,ret.end() , change_oo);
        return ret ;
    }

    std::set<std::string> reads ;

    static std::vector<ReadInfo> Add( const std::vector<ReadInfo> & c1
            , const std::vector<ReadInfo> & c2 , bool o1 , bool o2  ) 
    {
        std::vector<ReadInfo> ret ;
        if( o1 && o2 )
        {
            std::copy(c1.begin(),c1.end(),std::back_inserter(ret)) ;
            std::copy(c2.begin(),c2.end(),std::back_inserter(ret)) ;
        }
        else if ( o1 && !o2 )
        {
            auto c2p = CompleteReverse(c2);
            std::copy(c1.begin(),c1.end(),std::back_inserter(ret)) ;
            std::copy(c2p.begin(),c2p.end(),std::back_inserter(ret)) ;
        }
        else if ( !o1 && !o2 )
        {
            auto c1p = CompleteReverse(c1);
            auto c2p = CompleteReverse(c2);
            std::copy(c1p.begin(),c1p.end(),std::back_inserter(ret)) ;
            std::copy(c2p.begin(),c2p.end(),std::back_inserter(ret)) ;
        }
        else
        {
            auto c1p = CompleteReverse(c1);
            std::copy(c1p.begin(),c1p.end(),std::back_inserter(ret)) ;
            std::copy(c2.begin(),c2.end(),std::back_inserter(ret)) ;
        }
        return ret ;
    }

    static std::vector<ReadInfo> FilterCommon( const std::vector<ReadInfo> & base , 
            const std::set<std::string> & common )
    {
        auto check = [&]( const ReadInfo & info) -> bool 
        {
            return common.find( info.order) != common.end() ;
        };
        std::vector<ReadInfo> ret ;
        std::copy_if( base.begin() , base.end() , ret.begin() , check) ;
        return ret ;
    }
};

std::set<std::string> SetUnion( const std::set<std::string> & s1 ,
        const std::set<std::string> & s2 )
{
    std::set<std::string> dest1;
    std::set_union(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

std::set<std::string> SetDiff( const std::set<std::string> & s1 ,
        const std::set<std::string> & s2 )
{
    std::set<std::string> dest1;
    std::set_difference(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

std::set<std::string> SetIntersection( const std::set<std::string> & s1 ,
        const std::set<std::string> & s2 )
{
    std::set<std::string> dest1;
    std::set_intersection(s1.begin(), s1.end(),
            s2.begin(), s2.end(),
            std::inserter(dest1,dest1.begin()));
    return dest1 ;
}

struct AppConfig 
{
    std::map<std::string , std::set<std::string> > read_on_ont;

    std::map<std::string , MapperInfo> read2ont ;

    std::map<std::string , MapperInfo> read2con ;

    std::string r2ont_f ;
    std::string r2con_f ;

    void  LoadR2ONT() 
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(r2ont_f);
        if( in == NULL )
            FATAL( " failed to open r2ont for read !! exit ..." );
        auto parseline = [this]( const std::string & line ) 
        {
            MapperInfo tmp ;
            tmp.InitFromStr(line);
            read2ont[tmp.ref_id] = tmp ;
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete in ;
    }
    void  LoadR2CON() 
    {
        auto in = BGIQD::FILES::FileReaderFactory
            ::GenerateReaderFromFileName(r2con_f);
        if( in == NULL )
            FATAL( " failed to open r2con for read !! exit ..." );
        auto parseline = [this]( const std::string & line ) 
        {
            MapperInfo tmp ;
            tmp.InitFromStr(line);
            read2con[tmp.ref_id] = tmp ;
        };
        BGIQD::FILES::FileReaderFactory::EachLine(*in,parseline);
        delete in ;
    }

    void BuildR2ONTIndex()
    {
        for( const auto & pair : read2ont )
        {
            for( const auto & read_i : pair.second.reads)
            {
                read_on_ont[read_i].insert(pair.first) ;
            }
        }
    }

    BGIQD::stLFR::ScaffInfoHelper helper ;

    void LoadScaffInfos()
    {
        helper.LoadAllScaff(std::cin);
    }

    bool ParseAGap(BGIQD::stLFR::ContigDetail & c1
            ,const BGIQD::stLFR::ContigDetail & c2 )
    {
        // step 1 , find all common ONT reads
        const auto &  c1info = read2con.at(std::to_string(c1.contig_id)) ;
        const auto &  c2info = read2con.at(std::to_string(c2.contig_id)) ;
        const auto &  rc1 = c1info.reads ;
        const auto &  rc2 = c2info.reads ;

        std::set<std::string> ont_c1 ;
        std::set<std::string> ont_c2 ;

        std::map<std::string , std::set<std::string> > ont_share_reads;

        for( const auto & r1 : rc1 )
        {
            auto itr = read_on_ont.find(r1);
            if( itr != read_on_ont.end() )
            {
                const auto & onts = itr->second ;
                for( const auto & ont : onts  )
                {
                    ont_share_reads[ont].insert(r1);
                    ont_c1.insert( ont );
                }
            }
        }

        for( const auto & r1 : rc2 )
        {
            auto itr = read_on_ont.find(r1);
            if( itr != read_on_ont.end() )
            {
                const auto & onts = itr->second ;
                for( const auto & ont : onts  )
                {
                    ont_share_reads[ont].insert(r1);
                    ont_c2.insert( ont );
                }
            }
        }
        auto common_ont = SetUnion(ont_c1 , ont_c2 );
        // step 2 , filter good ONT reads
        if( common_ont.empty() )
            return false ;
        // step 3 , align c1 + c2 to ONT reads 
        std::map<std::string , TheSmithWaterman> align_cache ;
        for( const auto & ont : common_ont )
        {
            const auto & common = ont_share_reads.at(ont) ;
            const auto & ont_reads_base = read2ont.at(ont).read_infos ;
            auto ont_common = MapperInfo::FilterCommon( ont_reads_base , common );

            const auto & c1_reads_base = read2con.at(c1info.ref_id).read_infos;
            const auto & c2_reads_base = read2con.at(c2info.ref_id).read_infos;
            auto r1_common = MapperInfo::FilterCommon(c1_reads_base,common);
            auto r2_common = MapperInfo::FilterCommon(c2_reads_base,common);

            std::for_each(r1_common.begin() , r1_common.end() 
                    ,[c1](ReadInfo & i){ i.ref_id = std::to_string(c1.contig_id) ; });
            std::for_each(r2_common.begin() , r2_common.end() 
                    ,[c2](ReadInfo & i){ i.ref_id = std::to_string(c2.contig_id) ; });

            auto l2r = MapperInfo::Add( r1_common, r2_common 
                    , c1.orientation , c2.orientation);
            auto r2l = MapperInfo::CompleteReverse( l2r );

            BGIQD::stLFR::Schemes the_schemes ;
            the_schemes.match_orientation_score = 4 ;
            the_schemes.match_order_score = 2 ;
            the_schemes.delete_score = 0 ;
            the_schemes.insert_score = 0 ;
            the_schemes.match_orientation_score = 0 ;

            TheSmithWaterman tmp ;
            tmp.schemes = the_schemes ;
            tmp.ref = ont_common ;
            tmp.query = l2r ;
            tmp.CheckLoadedData() ;
            tmp.InitAfterDataLoaded();
            tmp.FillMutrix() ;

            TheSmithWaterman tmp1 ;
            tmp1.schemes = the_schemes ;
            tmp1.ref = ont_common ;
            tmp1.query = r2l ;
            tmp1.CheckLoadedData() ;
            tmp.InitAfterDataLoaded();
            tmp1.FillMutrix() ;

            if( tmp.max_value > tmp1.max_value )
            {
                align_cache[ont]=tmp;
            }
            else
            {
                align_cache[ont]=tmp1;
            }
        }
        // step 5 , choose the best ONT read 
        int max_value = 0;
        TheSmithWaterman max_align ;
        for(const auto & pair : align_cache )
        {
            auto path = pair.second.GetResult();
            auto path_ref = pair.second.AlignedElementsRef(path);
            auto path_query = pair.second.AlignedElementsRef(path);
            BGIQD::FREQ::Freq<std::string> freq ;
            std::for_each( path_query.begin() , path_query.end() ,
                    [&freq]( ReadInfo & i ) { freq.Touch(i.ref_id) ; } );
            if( freq.data.size() < 2 )
                continue ;
            if( pair.second.max_value > max_value )
            {
                max_value = pair.second.max_value;
                 max_align = pair.second ;
            }
        }
        align_cache.clear() ;
        if( max_value < 1 )
            return  false ;
        // step 6 , estimate the gap size.
        auto path       = max_align.GetResult();
        auto path_ref   = max_align.AlignedElementsRef(path);
        auto path_query = max_align.AlignedElementsRef(path);

        auto prev = std::adjacent_find(path_query.begin() 
                , path_query.end() 
                , [](const ReadInfo & p , const ReadInfo & n)
                        { return p.ref_id != n.ref_id ; } );

        if ( prev == path_query.end() )
            return false ;
        auto prev_from_query = *prev ;
        auto next_from_query = *( prev +1 ) ;
        auto prev_from_ref = *(path_ref.begin() + ( prev - path_query.begin() ) );
        auto next_from_ref = *(path_ref.begin() + ( prev - path_query.begin() +1 ) );
        // gap = ref_gap - c1_tail - c2_tail
        int ref_gap = std::abs( prev_from_ref.pos - next_from_ref.pos -1 );
        int c1_tail = 0 ;
        int c2_tail = 0 ;
        if( prev_from_ref.ref_id == std::to_string( c1.contig_id) )
        {
            c1_tail = c1.contig_len - prev_from_query.pos ;
            c2_tail = next_from_query.pos ;
        }
        else
        {
            c2_tail = c2.contig_len - prev_from_query.pos ;
            c1_tail = next_from_query.pos ;
        }
        int gap = ref_gap - c1_tail - c2_tail ;
        c1.gap_size = gap ;
        return true ;
    }

    void ParseAllGaps()
    {
        for( auto & pair : helper.all_scaff ) 
        {
            auto & a_scaff = pair.second.a_scaff ;
            for( int i = 0 ; i < (int) a_scaff.size() -1 ; i ++ )
            {
                ParseAGap(a_scaff.at(i) , a_scaff.at(i+1));
            }
        }
    }

    void PrintScaffInfos() 
    {
        helper.PrintAllScaff(std::cout);
    }
} config ;

int main(int argc , char **argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , r2ont, "read 2 ont data ");
        DEFINE_ARG_REQUIRED(std::string , r2con, "read 2 contig data ");
    END_PARSE_ARGS

    config.LoadR2ONT() ;
    config.BuildR2ONTIndex() ;
    config.LoadR2CON();
    config.LoadScaffInfos();
    config.ParseAllGaps() ;
    config.PrintScaffInfos();
    return 0;
}
