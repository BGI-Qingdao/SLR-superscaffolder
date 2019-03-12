#include <string>
#include <iostream>

#include "stLFR/CBB.h"
#include "common/freq/freq.h"

int main()
{
    BGIQD::FREQ::Freq<int> freq;
    std::string line;
    while(!std::getline(std::cin,line).eof())
    {
        BGIQD::stLFR::ContigBarcodeInfo tmp ;
        tmp.InitFromString(line);
        for( const auto x : tmp.barcodesOnPos)
        {
            freq.Touch( (int) x.second.size()) ;
        }
    }
    std::cout<<freq.ToString()<<std::endl;
    return 0;
}
