#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <tuple>

/******************************************************************************
 *
 * Basic functions 
 *
 *****************************************************************************/
inline size_t min(size_t i , size_t j) 
{
    return i < j ? i : j ;
}

inline size_t sum(const std::map<size_t , size_t> & hash)
{
    size_t ret = 0 ;
    for( auto pair : hash )
    {
        ret += pair.second;
    }
    return ret;
}

inline std::tuple<int,int,double> simularity(
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
    size_t all = sum(hash1) + sum(hash2) - both;
    double frac= double( both ) / double( all);
    return std::make_tuple(all,both,frac);
}

void median(std::map<size_t,double> & map , const std::map<size_t , size_t> index)
{
    for( auto itr = map.begin() ; itr != map.end() ; itr = std::next(itr) )
    {
        itr->second = itr->second / index.find(itr->first)->second ;
    }
}



/******************************************************************************
 *
 * Basic Marco and template
 *
 *****************************************************************************/

#define CHECK(x) _CHECK(x)
 #define _CHECK(x) \
     if(!(x)) {\
        std::cerr<<#x<<" failed at "<<__FILE__<<":"<<__LINE__<<std::endl;\
        exit(-1);\
     }

template<class T, class K ,class V>
void incr(T & map,K key, V v)
{
    auto itr = map.find(key);
    if( itr == map.end() )
    {
        map[key] = v ;
    }
    else
    {
        itr->second += v;
    }
}


template<class T>
void map2csv(const std::map<size_t , std::string> &heads
    , const T & map
    , size_t column 
    , size_t row
    , const std::string &file 
    , bool use_head = true 
    )
{
    std::ofstream ofs;
    ofs.open(file) ;
    if( ofs.bad() ) 
    {
        std::cerr<<"Failed to open "<<file<<std::endl;
        exit(-1);
    }
    else
    {
        std::cerr<<"Succ to open "<<file<<" for write"<<std::endl;
    }
    if( use_head )
    {
        for ( size_t j = 1 ; j<= row ; j++ )
        {
            ofs<<heads.at(j) <<",";
        }
        ofs<<std::endl;
    }
    for ( size_t c = 1 ; c <= column ; c++ )
    {
        for ( size_t r = 1 ; r<= row ; r++ )
        {
            ofs<<map.at(r).at(c)<<",";
        }
        ofs<<std::endl;
    }
    ofs.close();
}

size_t random_from_map(std::map<size_t,size_t> & map , size_t & max)
{
    size_t r = std::rand() % max ;
    size_t next = 0;
    for( auto & pair : map) 
    {
        next += pair.second;
        if(next > r)
        {
            pair.second -= 1;
            max -=1 ;
            if ( pair.second == 0 )
            {
                map.erase(pair.first) ;
            }
            return pair.first;
        }
    }
    CHECK(false);
    return -1;
}
void shuffle_map(std::map<size_t , std::map<size_t,size_t> > & map)
{
    std::srand(std::time(0));
    std::map<size_t, size_t> pool;
    // collect  all barcode 
    for( auto & pair1 : map )
    {
        for( auto pair2: pair1.second )
        {
            incr (pool ,pair2.first, pair2.second);
        }
    }
    size_t pool_size = sum( pool ) ;
    std::cerr<<"Total "<<pool_size<<" barcodes"<<std::endl;
    size_t count = 0 ;
    // shuffle barcode back.
    for ( auto & pair1 : map )
    {
        auto & a_bin = pair1.second;
        size_t bin_barcode = sum(a_bin);
        a_bin.clear();
        for( int i = 0 ; i < bin_barcode ; i++ )
        {   
            count++;
            if( count % 100== 0 )
            {
                std::cerr<<"shuffle done "<<count<<" barcodes"<<std::endl;
            }
            incr(a_bin, random_from_map(pool,pool_size) ,1);
        }
    }
    CHECK(sum(pool) == 0 );
}
/******************************************************************************
 *
 * Basic test  
 *
 *****************************************************************************/
void Test()
{
    std::cerr<<"Start basic test ..."<<std::endl;
    CHECK( min(1,2) == 1);
    std::map<size_t , size_t> hash1;
    hash1[1]=1;
    hash1[2]=1;
    hash1[3]=1;
    std::map<size_t , size_t> hash2;
    hash2[1]=1;
    hash2[2]=2;
    CHECK(sum(hash1) ==3 );
    CHECK(sum(hash2) ==3 );
    CHECK(simularity(hash1,hash2) ==  std::make_tuple(4,2,0.5f));
    std::cerr<<"End basic test ..."<<std::endl;

}

/******************************************************************************
 *
 * Logic.
 *
 *****************************************************************************/

int main(int argc , char ** argv)
{
    // Test first ...
    Test();   

    // Parse commands ...
    std::string input_file;
    std::string output_file;
    int step_max = 0;
    int limit= 0;
    int calc_step;
    int opt = 0;
    while( ( opt =getopt(argc,argv,"i:s:b:a:") ) !=EOF )
    {
        char data[1024];
        switch(opt)
        {
            case 'i':
                sscanf(optarg,"%s",data);
                input_file= std::string(data);
                break;
            case 'b':
                sscanf(optarg,"%d",&step_max);
                break;
            case 'a':
                sscanf(optarg,"%d",&limit);
                break;
            case 's':
                sscanf(optarg,"%d",&calc_step);
                break;
            default:
                break;
        }
    }

    const std::string Usage("Usage: bin_simularity -i ifile  -b bin_size -s calc_step  [ -a row_max ][-d]");
    if ( input_file == "" || step_max < 1 || calc_step < 1)
    {
        std::cerr<<Usage<<std::endl;
        exit(-1);
    }

    // Loading data ...
    std::ifstream ifs;
    ifs.open(input_file);
    if(ifs.bad())
    {
        std::cerr<<"Failed to open "<<input_file<<std::endl;
        exit(-1);
    }
    else
    {
        std::cerr<<"Succ to open "<<input_file<<" for read"<<std::endl;
    }
    int bin_max = 0;
    std::map<size_t, std::map<size_t,size_t> > bin_data;
    while( ! ifs.eof() )
    {
        int bin,num;
        char dot;
        ifs>>bin>>dot>>num;
        //std::cerr<<"DEBUG : "<<barcode<<":"<<num<<std::endl;
        if ( limit != 0 && bin>limit +step_max )
        {
            break;
        }
        if( bin > bin_max )
        {
            bin_max = bin ;
        }
        for(int i = 0 ; i < num ; i++ )
        {
            size_t key, value;
            ifs>>key>>dot>>value;
            bin_data[bin][key] = value;
            //std::cerr<<"DEBUG : "<<key<<":"<<value<<std::endl;
        }
    }
    ifs.close();

    std::map<size_t , size_t> total_index;
    std::map<size_t , std::map<size_t , double> > total_total;
    std::map<size_t , std::map<size_t , size_t> > detail_U;
    std::map<size_t , std::map<size_t , size_t> > detail_N;
    std::map<size_t , std::map<size_t , float> >  detail_S;
    if ( limit == 0 || limit > bin_max-step_max )
    {
        limit = bin_max-step_max;
    }
    // Calc now .,..
    for( int i = 0 ; i <= limit; i++ )
    {
        if( i % 1000 == 0 )
        {
            std::cerr<<"process "<<i<<" bin ..."<<std::endl;
        }
        auto itr1 = bin_data.find(i);
        if( itr1 == bin_data.end() )
        {
            std::cout<<0<<std::endl;
            continue;
        }
        int j = i + calc_step ;
        //for(int j = 1 ; j<=step_max ; j++ )
        {
            auto itr2 = bin_data.find(j) ;
            if( itr2 == bin_data.end() )
            {
                std::cout<<0<<std::endl;
                continue;
            }
            auto ret = simularity( itr1->second , itr2->second );
            std::cout<<std::get<2>(ret)<<std::endl;
        }
    }
    return 0;
}
