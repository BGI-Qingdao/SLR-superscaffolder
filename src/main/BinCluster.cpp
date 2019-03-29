#include "algorithm/incr_array/incr_array.h"
#include "algorithm/collection/collection.h"

#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/multithread/MultiThread.h"
#include "common/stl/mapHelper.h"
#include "common/log/log.h"
#include "common/error/Error.h"
#include "common/log/logfilter.h"
#include "common/files/file_reader.h"
#include "common/middle_valiad/MiddleValid.h"
#include "common/freq/freq.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/CBB.h"
#include "stLFR/TrunkGap.h"

#include <algorithm>
#include <iostream>
#include <string>

struct AppConfig
{
    typedef std::map< int , std::set<int> > BinIndexOnBarcode;

    typedef std::map<BGIQD::SOAP2::ContigId ,std::map< BGIQD::SOAP2::ContigId ,float> > ContigSims;

    BGIQD::stLFR::BarcodeOnBinArray barcodeOnBin ;

    BGIQD::stLFR::BinRelationArray  relations;

    ContigSims contigSims;

    BGIQD::stLFR::ContigRelationArray contig_relations;

    BinIndexOnBarcode binOnBarcode ;

    float thresold;


    BGIQD::SOAP2::FileNames fName;

    BGIQD::LOG::logger lger;

    bool calc_same_contig;


    enum WorkMode
    {
        Unknow = 0 ,
        Jaccard = 1 ,
        JoinBarcodeNum = 2
    } ;

    WorkMode work_mode ;

    void Init(const std::string & p , float t, bool cs)
    {
        fName.Init(p);
        thresold = t;
        calc_same_contig = cs ;
        barcodeOnBin.Init(1024);
        relations.Init(1024);
        contig_relations.Init(1024);
        BGIQD::LOG::logfilter::singleton().get("BinCluster",BGIQD::LOG::INFO,lger);
        lger<<BGIQD::LOG::lstart()<<"Init finsish ..."<<BGIQD::LOG::lend();
    }

    void LoadB2BArray( )
    {
        BGIQD::stLFR::LoadBarcodeOnBinArray( fName.BarcodeOnBin(middle_name) , barcodeOnBin);
    }

    void AllocRelationArray()
    {
        relations.init_n_element(barcodeOnBin.size());
    }
    std::set<size_t> del_indexs;
    void BuildBinIndexOnBarcode()
    {
        std::vector<int> size_cache;
        for(size_t i = 0 ; i <barcodeOnBin.size() ; i++)
        {
            auto & b2b = barcodeOnBin.at(i);
            for( auto j : b2b.collections )
            {
                int barcode = j.first ;
                binOnBarcode[barcode].insert( i);
            }
        }
    }

    void Calc1Bin2All( int i )
    {
        auto & bin = barcodeOnBin[i];
        auto & result = relations[i];
        result.binId = bin.binId ;
        result.contigId = bin.contigId ;
        result.start = bin.start ;
        result.end = bin.end;
        result.binIndex = i ;
        if( ! IsBinMiddleValid( bin.collections ) )
            return ;

        if( ! same_bin_only )
        {
            CalcSimForAllRelateBin( bin , result , i );
        }
        else
        {
            CalcSimForSameBin( bin , result , i );
        }
    }

    void CalcSimForAllRelateBin( const BGIQD::stLFR::BarcodeOnBin & bin 
            , BGIQD::stLFR::BinRelation & result
            , int i_index )
    {
            std::set<int> relates ;
            for(auto pair : bin.collections)
            {
                int barcode = pair.first ;
                if( binOnBarcode[barcode].size() > 10000 )
                {
                    continue ;
                }
                for( auto index : binOnBarcode[barcode] )
                {
                    if(! calc_same_contig && barcodeOnBin[index].contigId == bin.contigId )
                        continue ;
                    if ( index > i_index )
                        relates.insert(index);
                }
            }
            CalcSimWithRelates(bin,result,relates );
    }

    void CalcSimForSameBin( const BGIQD::stLFR::BarcodeOnBin & bin 
            , BGIQD::stLFR::BinRelation & result
            , int i_index )
    {
        std::set<int> relates ;
        // Only deal with same contig bin
        for(auto pair : bin.collections)
        {
            int barcode = pair.first ;
            for( auto index : binOnBarcode[barcode] )
            {
                if ( barcodeOnBin[index].contigId == bin.contigId && index > i_index )
                    relates.insert(index);
            }
        }
        CalcSimWithRelates(bin,result,relates );
    }

    void CalcSimWithRelates( const BGIQD::stLFR::BarcodeOnBin & bin 
            , BGIQD::stLFR::BinRelation & result 
            , const std::set<int> &  relates )
    {
        for(auto index : relates)
        {
            auto & other = barcodeOnBin[index];
            if( ! IsBinMiddleValid( bin.collections ) )
                continue;
            if( neib_only && ( !IsNeibValid ( bin , other ) ) )
                continue ;
            float sim = 0 ;
            if( work_mode == WorkMode::Jaccard )
                sim = BGIQD::stLFR::BarcodeCollection::Jaccard(bin.collections,other.collections);
            else if ( work_mode == WorkMode::JoinBarcodeNum )
                sim =  BGIQD::stLFR::BarcodeCollection::Intersection(bin.collections,other.collections).size();
            else
                assert(0);

            if( sim >= thresold )
            {
                auto & sinfo = result.sims[index];
                sinfo.binIndex = index ;
                sinfo.contigId = other.contigId;
                sinfo.binId = other.binId;
                sinfo.simularity= sim;
            }
        }
    }

    void RunAllJob( int thread)
    {
        BGIQD::MultiThread::MultiThread t_jobs;
        t_jobs.Start(thread);

        for(size_t i = 0 ; i <barcodeOnBin.size() ; i++)
        {
            t_jobs.AddJob([i,this](){ Calc1Bin2All(i);});
        }

        t_jobs.End();
        t_jobs.WaitingStop();
    }

    void CleanCollectionAndIndexMap()
    {
        barcodeOnBin.deep_clean();
        binOnBarcode.clear();
    }
    void BuildABBAResult()
    {
        for( size_t i = 0 ; i < relations.size() ; i++ )
        {
            auto & result = relations.at(i);
            for(auto pair : result.sims )
            {
                if( (int)pair.first <=(int) i )
                    continue ;
                auto & sinfo = pair.second ;
                auto & another = relations.at(sinfo.binIndex);
                auto & new_sinfo = another.sims[i];
                new_sinfo.binIndex = i ;
                new_sinfo.contigId = result.contigId;
                new_sinfo.binId = result.binId ;
                new_sinfo.simularity = sinfo.simularity ;
            }
        }
    }

    void PrintBinRelation()
    {
        BGIQD::stLFR::PrintBinRelationArray(fName.bin_cluster(middle_name) , relations);
    }

    void BuildContigRelation()
    {
        for( size_t i = 0 ; i < relations.size() ; i++ )
        {
            auto & result = relations.at(i);
            for(auto pair : result.sims )
            {
                auto & sinfo = pair.second ;
                BGIQD::STL::MapHelper<ContigSims::mapped_type>::UpdateAsBiggest(
                        contigSims[result.contigId]
                    ,   sinfo.contigId
                    ,   sinfo.simularity );
            }
        }
        for( auto pair : contigSims )
        {
            BGIQD::stLFR::ContigRelation data;
            data.contigId = pair.first ;
            for( auto & pair_2 : pair.second )
            {
                if( pair_2.first == pair.first )
                    continue ;
                BGIQD::stLFR::ContigSimularity sinfo;
                sinfo.contigId = pair_2.first ;
                sinfo.simularity = pair_2.second ;
                data.sims[pair_2.first] = sinfo;
            }
            contig_relations.push_back(data);
        }
    }

    void PrintContigRelation()
    {
        BGIQD::stLFR::PrintContigRelationArray(fName.cluster(middle_name) , contig_relations);
    }

    float del_fac ;
    bool check_freq_valid ;
    long biggest_freq ;
    long smallest_freq ;
    BGIQD::FREQ::Freq<int> bin_size_freq ;


    bool IsBinMiddleValid(
            const BGIQD::stLFR::BarcodeCollection & binc )
    {
        if( ! check_freq_valid )
            return true ;
        int freq = bin_size_freq.GetFreq(binc.keysize());
        if( freq < smallest_freq || freq > biggest_freq )
            return false ;
        return true ;
    }

    void BuildMiddleValid()
    {
        if( del_fac < 0.0000001f )
        {
            check_freq_valid = false ;
            biggest_freq = 0 ;
            smallest_freq = 0 ;
            return ;
        }
        check_freq_valid = true ;
        float valid = 1.0f - del_fac ;
        for( const auto & pair : barcodeOnBin )
        {
            bin_size_freq.Touch(pair.collections.keysize()) ;
        }
        std::tie(smallest_freq,biggest_freq) = 
            BGIQD::MIDDLE_VALID::MiddleValid(bin_size_freq.data,valid);
    }

    bool same_bin_only ; 

    std::string middle_name ;

    bool neib_only ;
    // contig_id <--> index , used to find neighbors.
    std::map<int , int > contig_indexs ;

    void LoadNeighbors()
    {
        if( ! neib_only )
            return ;

        auto in  = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(fName.mintreetrunklinear());
        if( in == NULL )
            FATAL(" failed to open xxx.mintree_trunk_linear for read!!! ");
        typedef BGIQD::stLFR::TrunkGap<int> GapInfo;
        std::map<int,std::vector<GapInfo>> gaps;
        BGIQD::stLFR::Load_MST_Trunk_Linear(*in, gaps);
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
        delete in ;
    }


    bool IsNeibValid(const BGIQD::stLFR::BarcodeOnBin bin 
            , const BGIQD::stLFR::BarcodeOnBin other)
    {
        if( ! neib_only )
            return true;
        if( contig_indexs.find( bin.contigId ) == contig_indexs.end()) 
            return false ;
        if( contig_indexs.find( other.contigId ) == contig_indexs.end()) 
            return false ;
        return std::abs( contig_indexs.at(bin.contigId) -
                contig_indexs.at(other.contigId) ) < 2;
    }

} config;

int main(int argc ,char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix. Input xxx.barcodeOnBin ; Output xxx.bin_cluster && xxx.cluster");
    DEFINE_ARG_OPTIONAL(std::string,  middle_name, "the middle name of output suffix " ,"");
    DEFINE_ARG_REQUIRED(float , threshold, "simularity threshold");
    DEFINE_ARG_OPTIONAL(int , thread, "thread num" ,"8");
    DEFINE_ARG_OPTIONAL(int , work_mode, "1 for Jaccard value , 2 for join_barcode_num " ,"1");
    DEFINE_ARG_OPTIONAL(float, del_fac, "del_fac or too small or too big bin set ." ,"0.0f");
    DEFINE_ARG_OPTIONAL(bool, pbc, "print bin cluster" ,"0");
    DEFINE_ARG_OPTIONAL(bool, bin_same_contig, "calc for bin on same contig ." ,"false");
    DEFINE_ARG_OPTIONAL(bool,  nb_only, "only calc bin for neighbor contig " ,"");
    DEFINE_ARG_OPTIONAL(bool, same_bin_only, "only calc for bin on same contig." ,"false");
    END_PARSE_ARGS

    config.del_fac = del_fac.to_float();
    config.neib_only = nb_only.to_bool() ;
    config.middle_name = middle_name.to_string() ;
    config.work_mode = static_cast<AppConfig::WorkMode>(work_mode.to_int());
    config.same_bin_only = same_bin_only.to_bool() ;
    config.Init(prefix.to_string() , threshold.to_float(), bin_same_contig.to_bool());
    BGIQD::LOG::timer t(config.lger,"BinCluster");

    config.LoadB2BArray() ;

    config.BuildBinIndexOnBarcode();

    config.AllocRelationArray();

    config.BuildMiddleValid() ;

    config.LoadNeighbors() ;

    config.RunAllJob(thread.to_int());

    config.CleanCollectionAndIndexMap() ;

    config.BuildABBAResult();

    if( pbc.to_bool() )
        config.PrintBinRelation();

    config.BuildContigRelation();

    config.PrintContigRelation();

    return 0;
}
