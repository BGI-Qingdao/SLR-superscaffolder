#ifndef __COMMON_FILES_FILE_READER_H__
#define __COMMON_FILES_FILE_READER_H__ 

/******************************************************************************
 *
 * This file provide user a module to open and read data from a file , without
 * special care about file type like txt or gzip .
 *
 *****************************************************************************/
#include <string>
#include <istream>

namespace BGIQD{
namespace FILES{

class FileReaderFactory {
    public:
        static std::istream* GenerateReaderFromFileName( const std::string & file_name );
        static void ResizeBuff(std::istream &, size_t size);
};

} //namespace FILES
} //namespace BGIQD
#endif //__COMMON_FILES_FILE_READER_H__
