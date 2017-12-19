#ifndef __COMMON_FILES_FILE_WRITER_H__
#define __COMMON_FILES_FILE_WRITER_H__

#include <ostream>
namespace BGIQD{
namespace FILES{

class FileWriterFactory{
    public:
        static std::ostream * GenerateWriterFromFileName( const std::string & file_name );

};

} //FILES
} //BGIQD

#endif //__COMMON_FILES_FILE_WRITER_H__
