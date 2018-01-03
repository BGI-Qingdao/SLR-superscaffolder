#include "contig_barcode.h"
#include "argsparser.h"
#include "file_reader.h"
#include "file_writer.h"
#include "sam_parser.h"
#include <cassert>

using namespace BGIQD;
using namespace BGIQD::JOB01;
using namespace BGIQD::ARGS;
using namespace BGIQD::LOG;
using namespace BGIQD::FILES;
using namespace BGIQD::SAM;


typedef std::map<int , int > contigTypes;

void loadContigTypes( const std::string & file , contigTypes & types )
{
    timer t(log1,"loadContigTypes");
    auto in = FileReaderFactory::GenerateReaderFromFileName(file);
    while( ! in->eof() )
    {
        std::string line;
        std::getline(*in,line);
        if( in->eof() )
            break;
        LineParser p(line);
        if(p.IsHead())
            continue;
        auto d0 = p.ParseAsMatchData();
        int read = std::stoi(d0.read_name);
        auto itr = types.find(read);
        if (d0.ref_name != "chr19" )
        {
            assert ( d0.ref_name == "*" );
            assert ( itr == types.end() );
            types[read] = 0 ; // means unmatched !
        }
        else
        {
            if( itr != types.end() )
            {
                assert( itr->second != 0);
                itr->second ++ ;
            }
            else
            {
                types[read] = 1;
            }
        }
    }
}

void printContigTypes( const std::string & file , const contigTypes & types)
{
    auto out = FileWriterFactory::GenerateWriterFromFileName(file);
    timer t(log1,"printContigTypes");
    for ( const auto & pair : types )
    {
        (*out)<<pair.first<<"\t"<<pair.second<<std::endl;
    }
    delete out;
}

int main(int argc ,char **argv)
{
    initLog("JOB05");

    START_PARSE_ARGS
    DEFINE_ARG(std::string , input , 'i');
    DEFINE_ARG(std::string , output, 'o');
    END_PARSE_ARGS

    contigTypes types;
    loadContigTypes(input.to_string() , types);
    printContigTypes(input.to_string() , types );
    return 0;
}
