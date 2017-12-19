#include "gzstream.h"
#include <fstream>
#include "file_writer.h"

namespace BGIQD{
namespace FILES{

    std::ostream * FileWriterFactory::GenerateWriterFromFileName( const std::string & file_name )
    {
        std::ostream * ret = NULL;
        if( file_name.rfind(".gz") == file_name.length() - 3 )
        {// gzip file
            ret = (new ogzstream(file_name.c_str()));
        }
        else
        { // treat as txt file now //TODO
            ret = (new std::ofstream(file_name.c_str()));
        }
        if (ret && ret->bad() )
        {
            delete ret ;
            ret = NULL;
        }
        return ret;
    }
}
}
