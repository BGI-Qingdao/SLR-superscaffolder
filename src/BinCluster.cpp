#include "contig_barcode.h"
#include "file_writer.h"
#include "argsparser.h"
#include <algorithm>
#include <iostream>
#include <common/multithread/MultiThread.h>
using namespace BGIQD;
using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;
using namespace BGIQD::LOG;
using namespace BGIQD::FILES;

//                length  contig  sd      bin
//typedef std::tuple<size_t,size_t,size_t, size_t> binIndex;

struct binIndex
{
    unsigned int contig;
    int bin;

    bool operator < ( const binIndex & a ) const 
    {
        return contig < a.contig ? true : bin<a.bin ;
    }
};
typedef std::vector<binIndex> binIndexs;

typedef std::vector<binIndex> binIndexs;

inline size_t min(size_t i , size_t j) 
{
    return i < j ? i : j ;
}

//
// not thread safe . do not edit sum_cache in multi-thread
//
static std::map<binIndex, size_t> sum_cache;
inline size_t sum(const  binIndex & index ,   const std::map<size_t , size_t> & hash)
{
    auto itr =  sum_cache.find( index ) ;
    if ( itr == sum_cache.end() )
    {
        size_t ret = 0 ;
        for( auto pair : hash )
        {
            ret += pair.second;
        }
        sum_cache.emplace( index , ret );
        //std::cerr<<"bin size "<<ret<<std::endl;
        return ret;
    }
    //std::cerr<<"bin size "<<itr->second<<std::endl;
    return itr->second;
}


inline std::tuple<int,int,double> simularity(
    const  binIndex & index1 ,
    const  binIndex & index2 ,
    const std::map<size_t , size_t> & hash1
  , const std::map<size_t , size_t> & hash2)
{
        size_t both=0;
        for( auto pair : hash1 )
        {
            auto itr2 = hash2.find(pair.first);
            if( itr2 != hash2.end() )
            {
                both += min( pair.second , itr2->second);
            }
        }
        //std::cerr<<" -- both"<<both<<std::endl;
        size_t all = sum(index1,hash1) + sum(index2,hash2) - both;
        //std::cerr<<" -- all "<<all<<std::endl;
        double frac= double( both ) / double(all);
        return std::make_tuple(all,both,frac);
}

void buildBinIndexs( const binBarcodeInfo & bbi , binIndexs & indexs)
{
    timer t(log1,"buildBinIndexs");
    for( const  auto & c : bbi )
    {
        for( const auto & b : c.second )
        {
            binIndex tmp ;
            tmp.contig = c.first;
            tmp.bin = b.first ;
            sum(tmp, b.second);
            indexs.push_back(tmp);//( c.first , b.first );
        }
    }
}

void updateMap( std::map<size_t,float> & data , size_t key ,float i )
{
    auto itr = data.find(key);
    if( itr == data.end() ) data[key] = i ;
    else if ( itr->second < i ) data[key] = i;
}

std::map<size_t,float>  cluster(const binIndex &seed , const binIndexs & allIndexs ,const binBarcodeInfo & bbi, float thresold )
{
    timer t(log1,"cluster");
    std::map<size_t,float>  rets;
    updateMap(rets, seed.contig , 1.0f);
    auto const & seedmap = bbi.at(seed.contig).at( seed.bin) ;// std::get<1>(seed)).at(std::get<3>(seed));

    for ( const auto & index : allIndexs)
    {
        auto const & curr = bbi.at(seed.contig).at(seed.bin);
        auto ret = simularity( index , seed , curr , seedmap );
        if( std::get<2>(ret) >= thresold )
        {
            updateMap(rets,index.contig,std::get<2>(ret));
        }
    }
    return rets;
}

void printClusterData(const std::string & file ,const  std::vector< std::map< size_t ,float > > &results)
{
    timer t(log1,"printClusterData");
    auto out = FileWriterFactory::GenerateWriterFromFileName(file);
    for( const auto r : results )
    {
        for( const auto p : r )
        {
            (*out)<<p.first<<":"<<p.second<<"\t";
        }
        (*out)<<std::endl;
    }
    delete out;
}
int main(int argc ,char **argv)
{
    initLog("BinCluster");
    timer t(log1,"binCluster");
    int t_num = 8;
    float thresold_f;
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , input , 'i',false,"barcodeOnBin");
    DEFINE_ARG_DETAIL(std::string , output, 'o',false,"output");
    DEFINE_ARG_DETAIL(int , thread, 't',true,"thread num [ default 8] ");
    DEFINE_ARG_DETAIL(float , thresold, 's',false,"simularity thresold");
    END_PARSE_ARGS
    if( thread.setted )
    {
       t_num = thread.to_int(); 
    }
    thresold_f = thresold.to_float();

    binBarcodeInfo bbi;
    loadBinBarcodeInfo(input.to_string() , bbi);
    binIndexs allIndexs;
    buildBinIndexs(bbi,allIndexs);


    std::mutex ret_mut;
    std::vector< std::map< size_t ,float > > results;
    int   index = 0;

    auto cluster_job = [&results , &allIndexs ,&bbi ,
         thresold_f,&ret_mut,&index ]
             ( const binIndex & seed )
    {
        auto ret = cluster(seed,allIndexs,bbi,thresold_f);
        {
            std::lock_guard<std::mutex> l(ret_mut);
            results.push_back(ret);
            index ++ ;
            if( index % 1000 == 0 )
            {
                log1<<BGIQD::LOG::lstart()<<"cluster "<<index<<" ... "<<BGIQD::LOG::lend();
            }
        }
    };
    BGIQD::MultiThread::MultiThread t_jobs;
    //TODO : make it thread safe
    t_jobs.Start(t_num);
    for( const auto & seed : allIndexs)
    {
        t_jobs.AddJob(std::bind(cluster_job,seed));
    }
    t_jobs.End();
    t_jobs.WaitingStop();

    printClusterData(output.to_string(),results);
    return 0;
}
