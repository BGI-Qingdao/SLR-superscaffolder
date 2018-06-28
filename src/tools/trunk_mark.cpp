#include <iostream>
#include <set>
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"

int main(int argc , char **argv)
{
    START_PARSE_ARGS 
    DEFINE_ARG_REQUIRED(std::string , seedLinear , " the seed linear file ");
    DEFINE_ARG_REQUIRED(std::string , trunkLinear , " the trunk linear file ");
    DEFINE_ARG_REQUIRED(std::string , outPrefix, " the out file prefix");
    END_PARSE_ARGS

    std::map<unsigned int , std::set<int> > seedPos;
    auto sf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(seedLinear.to_string());
    std::string line;
    int line_index = 1 ;
    while(!std::getline(*sf,line).eof())
    {
        auto items = BGIQD::STRING::split(line,"\t");
        //unsigned int seed = std::stoul(items[0]);
        //if( seedPos.find()
        seedPos[std::stoul(items[0])].insert( line_index) ;
        line_index ++ ;
    }
    delete sf ;

    auto tf = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(trunkLinear.to_string());
    std::ostream * out = NULL ;
    int trunk_index = 1 ;
    while(!std::getline(*sf,line).eof())
    {
        if(line[0] == '-')
        {
            if( out != NULL )
                delete out ;
            std::string of = outPrefix.to_string() + std::to_string(trunk_index);
            trunk_index ++ ;
            out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(of);
            continue;
        }
        unsigned int contig = std::stoul(line);
        (*out)<<contig;
        for( auto x : seedPos[contig])
            (*out)<<'\t'<<x;
        (*out)<<'\n';
    }
    delete tf ;
    return 0 ;
}
