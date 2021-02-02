#include "gzstream.h"
#include <fstream>
#include "file_writer.h"

/******************************************************************************
 *
 * This file provide user a module to open and write data from a file , without
 * special care about file type like txt or gzip .
 *
 *****************************************************************************/
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
        if (ret && ! ret->good() )
        {
            delete ret ;
            ret = NULL;
        }
        return ret;
    }

}
}
