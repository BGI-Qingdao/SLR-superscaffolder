#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/args/argsparser.h"


int main(int argc , char **argv)
{
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string, input,"input");
    DEFINE_ARG_REQUIRED(std::string, output , "output");
    END_PARSE_ARGS

    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(input.to_string());
    auto ou = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(output.to_string());
    BGIQD::FILES::FileReaderFactory::ResizeBuff(*in,1000000000);
    BGIQD::FILES::FileWriterFactory::ResizeBuff(*ou,1000000000);

    std::string line ;
    while( ! std::getline(*in,line).eof() )
    {
        *ou<<line<<std::endl;
    }


    delete in;
    delete ou;
    return 0;
}
