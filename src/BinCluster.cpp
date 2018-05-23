#include "algorithm/incr_array/incr_array.h"
#include "algorithm/collection/collection.h"

#include "common/files/file_writer.h"
#include "common/args/argsparser.h"
#include "common/multithread/MultiThread.h"

#include "soap2/soap2.h"
#include "soap2/fileName.h"

#include "stLFR/CBB.h"

#include <algorithm>
#include <iostream>
#include <string>

struct AppConfig
{
    typedef std::map< int , std::set<int> > BinIndexOnBarcode;

    BGIQD::stLFR::BarcodeOnBinArray barcodeOnBin ;
    BGIQD::INCRARRAY::IncrArray<BGIQD::stLFR::BinRelation>  relations;

    BinIndexOnBarcode binOnBarcode ;

    float thresold;

    BGIQD::SOAP2::FileNames fName;

    void Init(const std::string & p , float t)
    {
        fName.Init(p);
        thresold = t;
    }

    void LoadB2BArray( )
    {
        BGIQD::stLFR::LoadBarcodeOnBinArray( fName.BarcodeOnBin() , barcodeOnBin);
    }

    void AllocRelationArray()
    {
        relations.init_n_element(barcodeOnBin.size());
    }
    void BuildBinIndexOnBarcode()
    {
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

    void Calc1Bin2All(int i )
    {
        auto & bin = barcodeOnBin.at(i);
        auto & result = relations.at(i);
        result.binId = bin.binId ;
        result.contigId = bin.contigId ;
        result.binIndex = i ;
        std::set<int> relates ;
        for(auto pair : bin.collections)
        {
            int barcode = pair.first ;
            for( auto index : binOnBarcode[barcode] )
            {
                if ( index > i )
                    relates.insert(index);
            }
        }

        for(auto index : relates)
        {
            auto & other = barcodeOnBin.at(index);
            float sim = BGIQD::stLFR::BarcodeCollection::Jaccard(bin.collections,other.collections);
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

    void BuildABBAResult()
    {
        for( size_t i = 0 ; i < relations.size() ; i++ )
        {
            auto & binInfo = barcodeOnBin.at(i);
            auto & result = relations.at(i);
            for(auto pair : result.sims )
            {
                if( pair.first <=(int) i )
                    continue ;
                auto & sinfo = pair.second ;
                auto & another = relations.at(sinfo.binIndex);
                auto & new_sinfo = another.sims[i];
                new_sinfo.binIndex = i ;
                new_sinfo.contigId = binInfo.contigId;
                new_sinfo.binId = binInfo.binId ;
                new_sinfo.simularity = sinfo.simularity ;
            }
        }
    }

    void PrintBinRalation()
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(fName.bin_cluster());
        for( size_t i = 0 ; i < relations.size() ; i++ )
        {
            auto &result = relations.at(i);
            //TODO : new format
            (*out)<<result.contigId<<':'<<result.binId;
            for( auto pair : result.sims )
            {
                auto & sinfo = pair.second ;
                (*out)<<'\t'<<sinfo.contigId<<':'<<sinfo.simularity;
            }
            (*out)<<'\t'<<std::endl;
        }
        delete out;
    }
    /*
    void PrintContigRalation()
    {
        for( size_t i = 0 ; i < relations.size() ; i++ )
        {
            auto &result = relations.at(i);
            //TODO : new format
            (*out)<<result.contigId<<':'<<result.binId;
            for( auto pair : result.sims )
            {

            }
        }
    }*/

} config;

int main(int argc ,char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , prefix, "prefix");
    DEFINE_ARG_REQUIRED(float , thresold, "simularity thresold");
    DEFINE_ARG_OPTIONAL(int , thread, "thread num" ,"8");
    END_PARSE_ARGS

    config.Init(prefix.to_string() , thresold.to_float());

    config.LoadB2BArray() ;

    config.BuildBinIndexOnBarcode();

    config.AllocRelationArray();

    config.RunAllJob(thread.to_int());

    config.BuildABBAResult();

    config.PrintBinRalation();

    return 0;
}
