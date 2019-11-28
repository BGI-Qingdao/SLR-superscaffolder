#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/error/Error.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/freq/freq.h"
#include "soap2/fileName.h"
#include "stLFR/TrunkGap.h"
#include "stLFR/ScaffInfo.h"

struct Config {

    BGIQD::LOG::logger loger;
    struct GapExtra
    {
        bool prev_positive;
        bool next_positive;
        int  gap_len ;
        int  gap_orient ;  // 1 or 0 or -1 
        std::string ref ;
    };

    typedef BGIQD::stLFR::TrunkGap<GapExtra> GapInfo;
    std::map<int,std::vector<GapInfo>> gaps;

    void LoadTrunk( const std::string & file )
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        BGIQD::stLFR::Load_MST_Trunk_Linear(*in, gaps);
        delete in ;
        loger<<BGIQD::LOG::lstart() << "Load Trunk done "<<BGIQD::LOG::lend() ;
    }

    struct QuastInfo
    {
        int ref_start ;
        int ref_end ;
        int query_start ;
        int query_end ;
        std::string ref ;
        unsigned int contig ;
        bool orient ;
        int rank ;
        void InitFromString(const std::string & line )
        {
            std::string tmp1 ;
            std::istringstream ist(line);
            ist>>ref_start>>ref_end>>query_start>>query_end>>ref>>tmp1;
            std::string ctg ;
            for( auto i : tmp1 )
            {
                if( std::isdigit(i) )
                {
                    ctg += i;
                }
                else
                    break ;
            }
            contig = std::stoul(ctg);
            if( query_end < query_start )
                orient = false ;
            else 
                orient = true ;
        }
    };

    std::map<unsigned int ,QuastInfo > contigs ;

    void LoadSorted(const std::string & file )
    {
        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
        if( in == NULL )
            FATAL(" failed to open sorted_unique for read!!! ");
        std::string line ;
        unsigned int rank = 1 ;
        while( ! (*in).eof() )
        {
            std::getline((*in),line);
            if(line.empty() ) continue ;
            QuastInfo tmp ;
            tmp.InitFromString(line);
            tmp.rank = rank ++ ;
            contigs[tmp.contig] = tmp;
        }
        delete in;
        loger<<BGIQD::LOG::lstart() << "Load sorted unique done "<<BGIQD::LOG::lend() ;
    }

    void ProcessPair() {
        for( auto & pair : gaps ) {
            auto & a_scaff = pair.second ;
            bool gap_valid = true ;
            for( auto & gap : a_scaff ) {
                if ( contigs.find( gap.prev ) == contigs.end() ) {
                    gap.data.prev_positive = true ;
                    gap_valid = false ;
                } else {
                    gap.data.prev_positive = contigs.at(gap.prev).orient ;
                }
                if ( contigs.find( gap.next ) == contigs.end() ) {
                    gap.data.next_positive = true ;
                    gap_valid = false ;
                } else {
                    gap.data.next_positive = contigs.at(gap.next).orient ;
                }
                if( ! gap_valid ) {
                    gap.data.gap_orient = 0 ;
                    gap.data.gap_len = 1000 ;
                    gap.data.ref = "";
                    continue ;
                }
                const auto prev_contig = contigs.at(gap.prev) ;
                const auto next_contig = contigs.at(gap.next) ;
                if( prev_contig.ref != next_contig.ref ) {
                    gap.data.gap_orient = 0 ;
                    gap.data.gap_len = 1000 ;
                    gap.data.ref = "";
                    continue ;
                }
                gap.data.ref = prev_contig.ref;
                if ( prev_contig.ref_end > next_contig.ref_end ) {
                    gap.data.gap_len = prev_contig.ref_start - next_contig.ref_end -1 ;
                    gap.data.gap_orient = -1 ;
                } else {
                    gap.data.gap_len = next_contig.ref_start - prev_contig.ref_end -1 ;
                    gap.data.gap_orient = 1 ;
                }
                if( std::abs(prev_contig.rank - next_contig.rank) > 2 && gap.data.gap_len > 30000 )
                    gap.data.gap_len = 30000 ;
            }
        }
        loger<<BGIQD::LOG::lstart() << "ProcessPair done "<<BGIQD::LOG::lend() ;
    }

    BGIQD::stLFR::ScaffInfoHelper helper ;
    void GenScaffInfos() {
        int scaff_id = 1 ;
        for( const auto & pair : gaps ){
            const auto & a_order = pair.second ;
            auto & a_scaff_info = helper.all_scaff[scaff_id];
            a_scaff_info.scaff_id = scaff_id ++ ;
            int scaff_order  = 0 ;
            for( const auto & gap : a_order ) 
                scaff_order += gap.data.gap_orient;
            if( scaff_order > 0 ) {
                for( const auto & gap : a_order ) {
                    BGIQD::stLFR::ContigDetail tmp ;
                    tmp.contig_id = gap.prev ;
                    tmp.gap_size = gap.data.gap_len ;
                    tmp.orientation = gap.data.prev_positive ;
                    a_scaff_info.a_scaff.push_back(tmp);
                }
                const auto & gap = *(a_order.rbegin()); 
                BGIQD::stLFR::ContigDetail tmp ;
                tmp.contig_id = gap.next ;
                tmp.orientation = gap.data.next_positive ;
                tmp.gap_size = 0 ;
                a_scaff_info.a_scaff.push_back(tmp);
            } else {
                for( auto  itr = a_order.rbegin() ; itr != a_order.rend() ; itr++ ) {
                    const auto & gap = *itr ;
                    BGIQD::stLFR::ContigDetail tmp ;
                    tmp.contig_id = gap.next ;
                    tmp.gap_size = gap.data.gap_len ;
                    tmp.orientation = gap.data.next_positive ;
                    a_scaff_info.a_scaff.push_back(tmp);
                }
                const auto & gap = *(a_order.begin()); 
                BGIQD::stLFR::ContigDetail tmp ;
                tmp.contig_id = gap.prev ;
                tmp.orientation = gap.data.prev_positive ;
                tmp.gap_size = 0 ;
                a_scaff_info.a_scaff.push_back(tmp);
            }
        }
        loger<<BGIQD::LOG::lstart() << "GenScaffInfos done "<<BGIQD::LOG::lend() ;
    }

    void PrintScaffInfos(const std::string & file) {
        auto  out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
        if( out == NULL )
            FATAL(" failed to open scaff_infos for write!!! ");
        helper.PrintAllScaff(*out);
        delete out ;
    }
    void Init() {
        BGIQD::LOG::logfilter::singleton().get("Order_ideal"
                ,BGIQD::LOG::loglevel::INFO,loger);
    }
} config;

int main(int argc , char **argv )
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(std::string , linear , "xxx.mintreetrunklinear");
        DEFINE_ARG_REQUIRED(std::string , sorted_unique , "sorted_unique_contig.txt");
        DEFINE_ARG_REQUIRED(std::string , scaff_infos , "xxx.scaff_infos");
    END_PARSE_ARGS
    config.Init();
    config.LoadTrunk(linear.to_string());
    config.LoadSorted(sorted_unique.to_string());
    config.ProcessPair();
    config.GenScaffInfos();
    config.PrintScaffInfos(scaff_infos.to_string());
    return 0 ;
}
