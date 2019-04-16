#include "algorithm/distribution/distribution.h"

#include "common/args/argsparser.h"

#include <sstream>
#include <string>

typedef BGIQD::DISTRIBUTION::IntervalDistribution<int> Items;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
        DEFINE_ARG_REQUIRED(int, start, " start ");
        DEFINE_ARG_REQUIRED(int, end, " end ");
        DEFINE_ARG_REQUIRED(int, bin_size, " bin_size ");
    END_PARSE_ARGS;
    Items dist;
    dist.Init(bin_size.to_int() , start.to_int() , end.to_int());
    std::string line;

    while( ! std::getline(std::cin,line).eof() )
    {
        int key , value ;
        std::istringstream ist(line);
        ist>>key>> value; 
        dist.Count(key,value);
    }
    dist.CalcPercent();
    std::cout<<dist.ToString()<<std::endl;
    return 0 ;
}
