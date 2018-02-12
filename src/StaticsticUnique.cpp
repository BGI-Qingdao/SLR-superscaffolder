#include "contig_barcode.h"
#include "argsparser.h"
#include <iostream>

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main()
{
    initLog("StatictisUnique");
    std::string line ;
    std::map<int , std::tuple<int , float , int > > contigs;

    long total = 0;
    double cov_all = 0.0f ;

    while(!std::getline(std::cin,line).eof())
    {
        int contigId , length;
        float cov ;
        int tip ;
        sscanf(line.c_str() ,">%d length %d cvg_%f_tip_%d",&contigId , &length , &cov , &tip);
        if( length > 2*63+1 && cov > 10 )
        {
            total += length;
            cov_all += (length * cov);
        }
        contigs[contigId] = std::make_tuple( length , cov , tip );
    }
    float Ecov = cov_all/(float) total;
    float UniqueLow = 0.5f*Ecov;
    float UniqueHigh = 1.5f*Ecov ;
    std::cout<<"E(cov) = "<< Ecov<< " range ( "<<UniqueLow<<" , "<<UniqueHigh<<")"<<std::endl;

    for( const auto &  contig : contigs)
    {
        float cov = std::get<1>(contig.second);
        if(cov > UniqueLow && cov < UniqueHigh )
        {
            std::cout<<contig.first<<"\t"<<"unique"<<std::endl;
        }
    }

    return 0;
}
