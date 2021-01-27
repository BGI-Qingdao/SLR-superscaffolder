#include "stLFR/StringIdCache.h"
#include "utils/files/file_reader.h"
#include "utils/files/file_writer.h"
#include "utils/string/stringtools.h"
#include "utils/error/Error.h"

#include <cassert>

namespace BGIQD{
namespace stLFR{

    long StringIdCache::Id(const std::string & tag)
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
            if(pair.size() <2 )
                continue ;
            data.AssignTag(pair.at(0), std::stoll(pair.at(1)));
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

    bool IdStringCache::HasId( long id ) 
    {
        if( data.find(id) != data.end() )
            return true;
        else
            return false ;
    }
    std::string IdStringCache::Id(long id)
    {
        if( data.find(id) != data.end() )
            return data.at(id) ;
        else
            return std::to_string(id);
    }

    void IdStringCache::LoadStringIdCache( const std::string & file)
    {
        if( ! preload )
            return ;
        auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
        if( in == NULL )
            FATAL( " open \"string <--> id\" map  file for read failed !!! " );
        std::string line ;
        while(in &&!std::getline(*in,line).eof())
        {
            auto pair = BGIQD::STRING::split(line,"\t");
            data[std::stoi(pair.at(1))]=pair.at(0);
        }
        if( in )
            delete in ;

    }

}
}
