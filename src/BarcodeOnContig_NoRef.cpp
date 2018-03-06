#include "common/files/file_reader.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/string/stringtools.h"
#include <vector>
#include <map>

using namespace BGIQD::ARGS;
BGIQD::LOG::logger log;

struct ConfigBarcodeInfo
{
    int length;
    std::map<unsigned int , std::vector<unsigned int > > barcodesOnPos;

    std::string format(unsigned int id) const 
    {
        std::ostringstream ost;
        ost<<id<<'\t'<<length;
        for( const auto & i : barcodesOnPos)
        {
            ost<<i.first<<':'<<i.second[0];
            for( int m = 1 ; m < (int)i.second.size(); m++)
            {
                ost<<'|'<<i.second[m];
            }
            ost<<'\t';
        }
        return ost.str();
    }
};

typedef std::map<unsigned int, ConfigBarcodeInfo >  BarcodeOnContig;


void LoadSeeds(const std::string & file , BarcodeOnContig & boc )
{
    BGIQD::LOG::timer t(log,"LoadSeeds");
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file) ;
    std::string line ;
    while( in && !std::getline(*in, line).eof() )
    {
        auto items = BGIQD::STRING::split(line,"\t");
        boc[std::stoul(items[0])].length = std::stoul(items[1]);
    }
    delete in ;
}


int main(int argc , char ** argv)
{

    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , seed, 's',false, "seed.txt");
    END_PARSE_ARGS
    BGIQD::LOG::logfilter::singleton().get("BarcodeOnContig_NoRef",BGIQD::LOG::loglevel::INFO,log);
    BarcodeOnContig boc;
    LoadSeeds(seed.to_string() , boc);

    {
        std::string line;   
        while(!std::getline(std::cin,line).eof())
        {
            long readId , contigId, pos , barcode;
            char dir;
            std::istringstream ist(line);
            ist>>readId>>contigId>>pos>>dir>>barcode;
            if( boc.find( contigId ) == boc.end() )
            {
                continue;
            }
            boc[contigId].barcodesOnPos[pos].push_back(barcode);
        }
        for( const auto & i : boc )
        {
            std::cout<<i.second.format(i.first)<<std::endl;
        }
    }
    return 0;
}
