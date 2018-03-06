#include "contig_barcode.h"
#include "file_writer.h"
#include "argsparser.h"
#include <algorithm>
#include <iostream>

using namespace BGIQD;
using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;
using namespace BGIQD::LOG;
using namespace BGIQD::FILES;

//                length  contig  sd      bin
typedef std::tuple<size_t,size_t,size_t, size_t> binIndex;
typedef std::vector<binIndex> binIndexs;

inline size_t min(size_t i , size_t j) 
{
    return i < j ? i : j ;
}

inline size_t sum(   const std::map<size_t , size_t> & hash)
{
        size_t ret = 0 ;
        for( auto pair : hash )
        {
            ret += pair.second;
        }
        //std::cerr<<"bin size "<<ret<<std::endl;
        return ret;
}

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

static std::map< std::pair<binIndex,binIndex> , std::tuple<int,int,double> > sim_cache;

inline std::tuple<int,int,double> simularity(
    const  binIndex & index1 ,
    const  binIndex & index2 ,
    const std::map<size_t , size_t> & hash1
  , const std::map<size_t , size_t> & hash2)
{
    //auto itr = sim_cache.find( std::make_pair(index1,index2) );
    //if ( itr == sim_cache.end() )
    //{
        size_t both=0;
        //std::cerr<<"simularity : "<<std::get<1>(index1)<<":"<<std::get<3>(index1)
        //    <<"\t"<<std::get<1>(index2)<<":"<<std::get<3>(index2);
        for( auto pair : hash1 )
        {
            auto itr2 = hash2.find(pair.first);
            if( itr2 != hash2.end() )
            {
                both += min( pair.second , itr2->second);
                //std::cerr<<"\t"<<pair.first ;
            }
        }
        //std::cerr<<" -- both"<<both<<std::endl;
        size_t all = sum(index1,hash1) + sum(index2,hash2) - both;
        //std::cerr<<" -- all "<<all<<std::endl;
        double frac= double( both ) / double(all);
        //sim_cache.emplace( std::make_pair( index1 , index2) , std::make_tuple(all,both,frac));
        return std::make_tuple(all,both,frac);
    //}
    //return itr->second;
}

void buildBinIndexs( const binBarcodeInfo & bbi , binIndexs & indexs)
{
    timer t(log1,"buildBinIndexs");
    size_t all1 = 0 , all_2 = 0;
    for( const  auto & c : bbi )
    {
        all1 ++;
        all_2 += c.second.size() ;
    }
    int mid1 = all_2/all1;
    for( const  auto & c : bbi )
    {
        int mid = c.second.size() /2 ;
        for( const auto & b : c.second )
        {
            indexs.emplace_back(3000000- abs(mid1 - (int)sum(b.second)) ,c.first, abs(int(b.first) - mid) , b.first );
        }
    }
    std::sort(indexs.begin() , indexs.end());
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
    updateMap(rets, std::get<1>(seed) , 1.0f);
    auto const & seedmap = bbi.at(std::get<1>(seed)).at(std::get<3>(seed));
    
    for ( const auto & index : allIndexs)
    //for( const  auto & c : bbi )
    {
        auto const & curr = bbi.at(std::get<1>(index)).at(std::get<3>(index));
        if( index < seed ) 
        {
            auto ret = simularity( index , seed , curr , seedmap );
            if( std::get<2>(ret) >= thresold )
            {
                updateMap(rets,std::get<1>(index),std::get<2>(ret));
            }
        }
        else
        {
            auto ret = simularity( seed , index ,  seedmap , curr );
            if( std::get<2>(ret) >= thresold )
            {
                updateMap(rets,std::get<1>(index),std::get<2>(ret));
            }
        }
    }
    for( const auto p : rets )
    {
        std::cout<<p.first<<":"<<p.second<<"\t";
    }
    std::cout<<std::endl;
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

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , input , 'i',false,"barcodeOnBin");
    DEFINE_ARG_DETAIL(std::string , output, 'o',false,"output");
    DEFINE_ARG_DETAIL(float , thresold, 's',false,"simularity thresold");
    END_PARSE_ARGS

    binBarcodeInfo bbi;
    loadBinBarcodeInfo(input.to_string() , bbi);
    binIndexs allIndexs;
    buildBinIndexs(bbi,allIndexs);


    std::vector< std::map< size_t ,float > > results;
    for( const auto & seed : allIndexs) 
    {
        results.push_back(cluster(seed,allIndexs, bbi, thresold.to_float()));
    }
    printClusterData(output.to_string(),results);
    return 0;
}
