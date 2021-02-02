#include "file_reader.h"
#include "gzstream.h"
#include <fstream>
#include <stdlib.h>
namespace BGIQD{
namespace FILES {

    std::istream * FileReaderFactory::GenerateReaderFromFileName(const std::string & file_name )
    {
        std::istream * ret = NULL;
        if( file_name.rfind(".gz") == file_name.length() - 3 )
        {// gzip file
            ret = (new igzstream(file_name.c_str()));
        }
        else
        { // treat as txt file now //TODO
            ret = (new std::ifstream(file_name.c_str()));
        }
        if (ret && ! ret->good() )
        {
            delete ret ;
            ret = NULL;
        }
        return ret;
    }

    void FileReaderFactory::EachLine(std::istream & ist , ParseLine parse)
    {
        std::string line;
        while( ! std::getline(ist,line).eof() )
        {
            parse(line);
        }
    }


}
}
