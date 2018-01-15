#include "common/files/file_reader.h"
#include "common/args/argsparser.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include <iostream>

using namespace BGIQD::ARGS;
using namespace BGIQD::LOG;
using namespace BGIQD::FILES;


int main(int argc ,char **argv)
{
    logger loger;
    logfilter::singleton().get("CHECK MAP OUTPUT", loglevel::INFO , loger);
    START_PARSE_ARGS
    DEFINE_ARG(std::string , input , 'i');
    END_PARSE_ARGS
    unsigned num_ctg = 406918384;
    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(input.to_string());
    long total =0 , err = 0 ;
    while(!in->eof() )
    {
        long long a ;
        unsigned int b ;
        int c;
        std::string line ;
        if((std::getline(*in,line)).eof())
            break;
        sscanf(line.c_str(),"%lld %u %d",&a,&b,&c);
        if ( b > num_ctg )
        {
            err ++;
        }
        total ++;
        if( total % 1000000000 )
            loger<<lstart()<<total<<" data pass ... "<<lend();
    }
    loger<<lstart()<<"---- summary ---- "<<lend();
    loger<<lstart()<<err <<" err in "<<total<<" data "<<lend();
    return 0;
}
