#include "common/args/argsparser.h"
#include "common/freq/freq.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"

#include "algorithm/linear_fitting/Minimum_multiplication.h"
#include "algorithm/collection/collection.h"

#include <iostream>
#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <tuple>

struct AppConfig
{
    typedef BGIQD::Collection::Collection<int> Bin;

    typedef BGIQD::LINEARFITTING::Item<int,float> LSItem;

    std::map<int , std::vector<float>> xys;

    std::vector<LSItem> xydata;

    int limit ;

    int step_max ;

    std::map<int,Bin> bin_data;

    int bin_max ;

    void LoadBarcodeOnRef()
    {
        while( ! std::cin.eof() )
        {
            int bin,num;
            char dot;
            std::cin>>bin>>dot>>num;
            //std::cerr<<"DEBUG : "<<barcode<<":"<<num<<std::endl;
            if ( limit != 0 && bin>limit +step_max )
            {
                break;
            }
            if( bin > bin_max )
                bin_max = bin ;
            for(int i = 0 ; i < num ; i++ )
            {
                size_t key, value;
                std::cin>>key>>dot>>value;
                bin_data[bin].IncreaseElement(key, value);
            }
        }
    }

    void CalcAll()
    {
        for( int i = 1 ; i <= bin_max - step_max ; i ++ )
        {
            for( int j = i + 1 ; j <= bin_max ; j++ )
            {
                xys[ j - i ].push_back( Bin::Jaccard(bin_data[i],bin_data[j]) );
            }
        }

        for( auto & x : xys )
        {
            float total = 0 ;
            for( float x : x.second )
            {
                total += x;
            }
            xydata.push_back( LSItem { x.first , total / x.second.size() } );
        }

    }
    void PrintLine()
    {
        auto ret = BGIQD::LINEARFITTING::lineFit(xydata);
        std::cout<<" a= "<<ret.a<<std::endl;
        std::cout<<" b= "<<ret.b<<std::endl;
        std::cout<<" c= "<<ret.c<<std::endl;
    }

} config ;

/******************************************************************************
 *
 * Logic.
 *
 *****************************************************************************/

int main(int argc , char ** argv)
{

    START_PARSE_ARGS 
        DEFINE_ARG_REQUIRED(int, limit, " limit ");
    DEFINE_ARG_REQUIRED(int, step, " step ");
    END_PARSE_ARGS;

    config.limit = limit.to_int();
    config.step_max = step.to_int();

    config.LoadBarcodeOnRef();
    config.CalcAll();
    config.PrintLine();

    return 0;
}
