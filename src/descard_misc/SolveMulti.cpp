#include <iostream>
#include <map>
#include <tuple>
#include <algorithm>
#include <vector>
#include <set>
#include "common/string/stringtools.h"
#include "common/args/argsparser.h"

typedef std::tuple<int,float,std::string> sort_t1;


std::map<std::string,std::vector<int> > map1;
std::map<std::string,std::vector<float> > map2;
std::map<std::string, std::map<int, sort_t1> > maps;
std::map<std::string , std::vector<std::tuple<float,std::string , int , float>>> map_value;
std::vector<int> len_map;
std::vector<float> sim_map;
std::vector<int> index_map;
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
        DEFINE_ARG_REQUIRED(int , type, "1 for len , 2 for simularity \
                                ,3 for split multi and base \
                                ,4 for solve multi");
        DEFINE_ARG_REQUIRED(int , len_max, "max_len for first close");
        DEFINE_ARG_REQUIRED(int , len_min, "min len for second close");
        DEFINE_ARG_REQUIRED(float, sim, "min simularity");
        DEFINE_ARG_REQUIRED(int, index, "index");
    END_PARSE_ARGS

    std::string line;
    int index1 = 0;
    while(!std::getline(std::cin,line).eof())
    {
        auto data = BGIQD::STRING::split(line,"\t");
        std::string key = data[0];
        std::string dir= data[1];
        for( size_t i = 2 ; i<data.size() ;i++)
        {
            int len ; float f ; std::string to ;
            auto ret = parse_string(data[i]);
            std::tie(len,f,to) = ret;
            map1[key+dir].push_back(std::get<0>(ret));
            map2[key+dir].push_back(std::get<1>(ret));
            maps[key+dir][std::get<0>(ret)] = ret;
            map_value[key+dir].push_back(std::make_tuple(len/f,to,len,f));
            len_map.push_back(len);
            sim_map.push_back(f);
            index1++;
        }
        index_map.push_back(index1);
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
    else if(type.to_int() == 3 )
    {
        int j = 0 ;
        for( size_t i = 0 ; i < len_map.size() ; i++ )
        {

            if ( (j == 0 && index_map[j]==1 )|| index_map[j] - index_map[j-1] == 1 ) 
            {
                std::cerr<<len_map[i]<<'\t'<<sim_map[i]<<std::endl;
            }
            else
            {
                std::cout<<len_map[i]<<'\t'<<sim_map[i]<<std::endl;
            }
            if( i + 1 >= size_t(index_map[j]) )
            {
                //std::cout<<"--"<<'\t'<<"--"<<std::endl;
                j++;
            }
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
        for( auto & i : map_value)
        {
            std::sort(i.second.begin() ,i.second.end());
        }
        std::map<int,std::set<int>> double_check;
        for( auto & i : map_value )
        {
            //if( i.second[0] <= len_max.to_int() && i.second[1] >= len_min.to_int() )
            {
                //auto & item  = maps[i.first][i.second[index.to_int()]];
                if( std::get<3>(i.second[0]) > sim.to_float() && std::get<0>(i.second[0]) * index.to_int() < std::get<0>(i.second[1]) )
                {
                    double_check[std::stoi(i.first)].insert(std::stoi(std::get<1>(i.second[0])));
                }
            }
        }
        std::set< std::pair<int,int> > dup;
        for( auto & i : double_check )
        {
            for( auto j : i.second )
            {
                if( dup.find( std::make_pair( i.first , j ) ) != dup.end() ) 
                        continue ;
                if(double_check.find(j) == double_check.end() || double_check[j].find(i.first) != double_check[j].end())
                {
                    std::cout<<i.first<<"\t"<<j<<std::endl;
                    dup.insert( std::make_pair( j , i.first ) );
                }
            }
        }
    }

    return 0;
}
