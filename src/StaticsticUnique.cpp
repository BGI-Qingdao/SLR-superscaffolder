#include "contig_barcode.h"
#include "argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include <iostream>

using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;


int main(int argc , char ** argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(int , kvalue, 'K',false,"K value");
    DEFINE_ARG_DETAIL(int , min, 'm',false,"min length of contig");
    END_PARSE_ARGS

    BGIQD::LOG::logger log;
    BGIQD::LOG::logfilter::singleton().get("StatictisUnique",BGIQD::LOG::loglevel::INFO,log);

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
        if( length > (2*kvalue.to_int())+1 && cov > 10 )
        {
            total += length;
            cov_all += (length * cov);
        }
        contigs[contigId] = std::make_tuple( length , cov , tip );
    }
    float Ecov = cov_all/(float) total;
    float UniqueLow = 0.5f*Ecov;
    float UniqueHigh = 1.5f*Ecov ;

    log<<BGIQD::LOG::lstart()<<"E(cov) = "<< Ecov<< " range ( "<<UniqueLow<<" , "<<UniqueHigh<<")"<< BGIQD::LOG::lend();
    long index = 0 ;
    long index1 = 0 ;
    for( const auto &  contig : contigs)
    {
        index1 ++ ;
        float cov = std::get<1>(contig.second);
        int len = std::get<0>(contig.second);
        if(cov > UniqueLow && cov < UniqueHigh && len > min.to_int() - kvalue.to_int())
        {
            index ++ ;
            std::cout<<contig.first<<'\t'<<len<<std::endl;
        }
    }

    log<<BGIQD::LOG::lstart()<<index << " seed in "<<index1<<" contigs"<< BGIQD::LOG::lend();

    return 0;
}
