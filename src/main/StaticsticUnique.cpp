#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "soap2/soap2.h"
#include <iostream>
#include <set>
#include <map>

struct AppConfig
{
    struct ContigInfo
    {
        BGIQD::SOAP2::ContigId id ;
        int length;
        float cov ;
        int  tip ;
    };

    int K;
    int min;
    long total ;
    float cov_all ;
    std::map<BGIQD::SOAP2::ContigId , ContigInfo> contigInfo;
    std::set<unsigned int> pass ;
    BGIQD::LOG::logger log;
    float Ecov ;
    float UniqueLow ;
    float UniqueHigh  ;

    void Init(int k, int m)
    {
        K = k ;
        min = m;
        total = 0;
        cov_all = 0 ;
        Ecov = 0;
        UniqueLow = 0 ;
        UniqueHigh = 0 ;
        contigInfo.clear();
        BGIQD::LOG::logfilter::singleton().get("StatictisUnique",BGIQD::LOG::loglevel::INFO,log);
    }

    void LoadContigInfio(std::istream & ist)
    {
        std::string line ;
        unsigned int prev = -1 ;
        while(!std::getline(ist,line).eof())
        {
            ContigInfo data;
            sscanf(line.c_str() ,">%u length %d cvg_%f_tip_%d",&data.id , &data.length , &data.cov , &data.tip);
            if( prev != (unsigned int )-1 && prev == data.id -1 )
                pass.insert(prev);
            prev = data.id ;
            if( data.length > (2*K)+1 && data.cov > 10 )
            {
                total += data.length;
                cov_all += (data.length * data.cov);
            }
            contigInfo[data.id] = data;
        }
    }

    void CalcCov()
    {
        Ecov = cov_all/(float) total;
        UniqueLow = 0.5f*Ecov;
        UniqueHigh = 1.5f*Ecov ;
        log<<BGIQD::LOG::lstart()<<"E(cov) = "<< Ecov<< " range ( "<<UniqueLow<<" , "<<UniqueHigh<<")"<< BGIQD::LOG::lend();
    }

    void PrintUniqueInfo(std::ostream & ost)
    {
        long index = 0 ;
        long index1 = 0 ;
        for( const auto &i : contigInfo)
        {
            auto & contig = i.second ;
            index1 ++ ;
            if( pass.find(contig.id) == pass.end() && contig.cov > UniqueLow && contig.cov < UniqueHigh && contig.length >=min )
            {
                index ++ ;
                ost<<contig.id<<'\t'<<contig.length<<std::endl;
            }
        }
        log<<BGIQD::LOG::lstart()<<index << " seed in "<<index1<<" contigs"<< BGIQD::LOG::lend();
    }

}config;

int main(int argc , char ** argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(int ,kvalue,"K value");
    DEFINE_ARG_REQUIRED(int ,min ," bin size ");
    END_PARSE_ARGS

    config.Init(kvalue.to_int(),min.to_int());
    BGIQD::LOG::timer t(config.log,"StaticsticUnique");
    config.LoadContigInfio(std::cin);
    config.CalcCov();
    config.PrintUniqueInfo(std::cout);

    return 0;
}
