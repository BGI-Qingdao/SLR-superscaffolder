#include "common/freq/freq.h"
#include "common/string/stringtools.h"

#include <iostream>
#include <stdio.h>
#include <cassert>

int main(int argc ,char ** argv)
{
    if( argc != 2 )
    {
        printf("Usage : %s  seed_len <seedLinear.txt\n",argv[0]);
        return 0;
    }
    int seed_len = std::atoi(argv[1]);
    std::cout<<"seed_len "<<seed_len<<std::endl;
    BGIQD::FREQ::Freq<int> freq;
    std::string line;
    int prev = -1 ;
    while( ! std::getline(std::cin,line).eof() )
    {
        auto item = BGIQD::STRING::split(line,"\t");
        assert(item.size() >= 4 );
        int len = std::stoi(item[1]);
        if( len < seed_len )
        {
            continue;
        }
        int s = std::stoi(item[2]);
        int e = std::stoi(item[3]);
        if( prev != -1 )
            freq.Touch((s - prev )/1000);
        prev = e ;
    }
    std::cout<<freq.ToString()<<std::endl;
    return 0 ;
}
