#include "algorithm/incr_array/incr_array.h"
#include "algorithm/collection/collection.h"

#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/multithread/MultiThread.h"
#include "common/stl/mapHelper.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/CBB.h"

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

    bool del;

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
        BGIQD::stLFR::LoadBarcodeOnBinArray( fName.BarcodeOnBin() , barcodeOnBin);
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
            if( del)
            {
                size_cache.push_back(b2b.collections.size());
            }
            for( auto j : b2b.collections )
            {
                int barcode = j.first ;
                binOnBarcode[barcode].insert( i);
            }
        }
        if( del )
        {

        }
    }

    void Calc1Bin2All(int i )
    {
        auto & bin = barcodeOnBin[i];
        auto & result = relations[i];
        result.binId = bin.binId ;
        result.contigId = bin.contigId ;
        result.start = bin.start ;
        result.end = bin.end;
        result.binIndex = i ;
        std::set<int> relates ;
        if( ! same_bin_only )
        {
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
                    if ( index > i )
                        relates.insert(index);
                }
            }
        }
        else
        {
            // Only deal with same contig bin
            for(auto pair : bin.collections)
            {
                int barcode = pair.first ;
                for( auto index : binOnBarcode[barcode] )
                {
                    if ( barcodeOnBin[index].contigId == bin.contigId && index > i )
                        relates.insert(index);
                }
            }
        }
        for(auto index : relates)
        {
            auto & other = barcodeOnBin[index];
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
        BGIQD::stLFR::PrintBinRelationArray(fName.bin_cluster() , relations);
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
        BGIQD::stLFR::PrintContigRelationArray(fName.cluster() , contig_relations);
    }

    bool same_bin_only ; 

} config;

int main(int argc ,char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix. Input xxx.barcodeOnBin ; Output xxx.bin_cluster && xxx.cluster");
    DEFINE_ARG_REQUIRED(float , threshold, "simularity threshold");
    DEFINE_ARG_OPTIONAL(int , thread, "thread num" ,"8");
    DEFINE_ARG_OPTIONAL(int , work_mode, "1 for Jaccard value , 2 for join_barcode_num " ,"1");
    DEFINE_ARG_OPTIONAL(bool, pbc, "print bin cluster" ,"0");
    DEFINE_ARG_OPTIONAL(bool, bin_same_contig, "calc for bin on same contig ." ,"false");
    DEFINE_ARG_OPTIONAL(bool, same_bin_only, "calc only for bin on same contig ." ,"false");
    //DEFINE_ARG_OPTIONAL(bool, del, "calc for bin on same contig ." ,"false");
    END_PARSE_ARGS

    config.work_mode = static_cast<AppConfig::WorkMode>(work_mode.to_int());
    config.same_bin_only = same_bin_only.to_bool() ;
    config.Init(prefix.to_string() , threshold.to_float(), bin_same_contig.to_bool());
    config.del = false;// del.to_bool();
    BGIQD::LOG::timer t(config.lger,"BinCluster");

    config.LoadB2BArray() ;

    config.BuildBinIndexOnBarcode();

    config.AllocRelationArray();

    config.RunAllJob(thread.to_int());

    config.CleanCollectionAndIndexMap() ;

    config.BuildABBAResult();

    if( pbc.to_bool() )
        config.PrintBinRelation();

    config.BuildContigRelation();

    config.PrintContigRelation();

    return 0;
}
