#include "stLFR/StringIdCache.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

namespace BGIQD{
namespace stLFR{

    int StringIdCache::Id(const std::string & tag)
    {
        if( tag == "0_0_0")
            return 0;
        if( preload )
            return data.Id(tag);
        else
        {
           return data.AddTag(tag); 
        }
    }

    void StringIdCache::Load( const std::string & file ) 
    {
        if( ! preload )
            return ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
        if( in == NULL )
            FATAL( " open (barcodeList) file for read failed !!! " );
        std::string line ;
        while(in &&!std::getline(*in,line).eof())
        {
            auto pair = BGIQD::STRING::split(line,"\t");
            data.AssignTag(pair.at(0), std::stoi(pair.at(1)));
        }
        if( in )
            delete in ;
    }
    void StringIdCache::Print(const std::string & file)
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
        if( out == NULL )
            FATAL( " open (barcodeList) file for write failed !!!");
        (*out)<<"0_0_0"<<'\t'<<0<<'\n';
        data.Print(*out);
        delete out;
    }
}
}
