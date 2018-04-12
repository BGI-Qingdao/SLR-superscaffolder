#include <iostream>
#include <map>
#include <tuple>
#include <algorithm>
#include <vector>
#include "common/string/stringtools.h"
#include "common/args/argsparser.h"

typedef std::tuple<int,float,std::string> sort_t1;


std::map<std::string,std::vector<int> > map1;
std::map<std::string,std::vector<float> > map2;
std::map<std::string, std::map<int, sort_t1> > maps;

sort_t1 parse_string(const std::string &s)
{
    auto items = BGIQD::STRING::split(s,":");
    return std::make_tuple(std::stoi(items[1]) , std::stof(items[2]),items[0]);
}

template<class T>
void printv(const std::vector<T> & t)
{
    for( const auto i : t)
    {
        std::cout<<i<<"\t";
    }
    std::cout<<std::endl;
}

int main(int argc , char ** argv)
{

    START_PARSE_ARGS
        DEFINE_ARG_DETAIL(int , type, 'o',false,"1 for len , 2 for simularity");
        DEFINE_ARG_DETAIL(int , len_max, 'a',false,"max_len for first close");
        DEFINE_ARG_DETAIL(int , len_min, 'i',false,"min len for second close");
        DEFINE_ARG_DETAIL(float, sim, 's',false,"min simularity");
        DEFINE_ARG_DETAIL(int, index, 'd',false,"index");
    END_PARSE_ARGS


    std::string line;
    while(!std::getline(std::cin,line).eof())
    {
        auto data = BGIQD::STRING::split(line,"\t");
        std::string key = data[0];
        std::string dir= data[1];
        for( size_t i = 2 ; i<data.size() ;i++)
        {
            auto ret = parse_string(data[i]);
            map1[key+dir].push_back(std::get<0>(ret));
            map2[key+dir].push_back(std::get<1>(ret));
            maps[key+dir][std::get<0>(ret)] = ret;
        }
    }

    if( type.to_int() == 1 )
    {
        for( auto & i : map1)
        {
            std::sort(i.second.begin() ,i.second.end());
            printv(i.second);
        }
    }
    else if( type.to_int() == 2 )
    {
        for( auto & i : map2)
        {
            std::sort(i.second.rbegin() ,i.second.rend());
            printv(i.second);
        }
    }
    else
    {
        for( auto & i : map1)
        {
            std::sort(i.second.begin() ,i.second.end());
        }
        for( auto & i : map2)
        {
            std::sort(i.second.rbegin() ,i.second.rend());
        }
        for( auto & i : map1 )
        {
            //if( i.second[0] <= len_max.to_int() && i.second[1] >= len_min.to_int() )
            {
                auto & item  = maps.at(i.first)[i.second[index.to_int()]];
                if( std::get<1>(item) >= sim.to_float() )
                {
                    std::cout<<std::stoi(i.first)<<"\t"<<std::get<2>(item)<<std::endl;
                }
            }
        }
    }

    return 0;
}
