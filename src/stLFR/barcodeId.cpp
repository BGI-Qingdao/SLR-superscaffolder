#include "stLFR/barcodeId.h"
#include "common/files/file_reader.h"
#include "common/files/file_writer.h"
#include "common/string/stringtools.h"
#include "common/error/Error.h"

namespace BGIQD{
namespace stLFR{

    BarcodeId BarcodeId::the_one ;
    int BarcodeId::AddTag(const std::string & tag)
    {
        auto itr = m_tag2num.find(tag);
        if( itr != m_tag2num.end() )
            return itr->second;
        while(1)
        {
            if( m_num2tag.find(curr) == m_num2tag.end() )
                break;
            curr ++ ;
        }
        m_num2tag[curr]=tag;
        m_tag2num[tag] = curr ;
        curr ++ ;
        return curr - 1;
    }

    int BarcodeId::Id(const std::string & tag ) const 
    {
        auto itr = m_tag2num.find(tag);
        if( itr != m_tag2num.end() )
            return itr->second;
        return -1;
    }

    bool BarcodeId::AssignTag( const std::string &tag, int number )
    {
        auto itr = m_tag2num.find(tag);
        if( itr != m_tag2num.end() )
            return itr->second == number ;
        auto itr1 = m_num2tag.find(number);
        if( itr1 != m_num2tag.end() )
            return false ;
        m_num2tag[number]=tag;
        m_tag2num[tag] = number;
        return true ;
    }

    void BarcodeId::Print( std::ostream & out ) const 
    {
        for( const auto & i : m_tag2num )
        {
            out<<i.first<<'\t'<<i.second<<'\n';
        }
    }

    /*********************************************************************/
    bool BarcodeIdHelper::preload = false ;
    int BarcodeIdHelper::Id(const std::string & tag)
    {
        if( tag == "0_0_0")
            return 0;
        if( preload )
            return BarcodeId::Singleton().Id(tag);
        else
        {
           return BarcodeId::Singleton().AddTag(tag); 
        }
    }

    void BarcodeIdHelper::Load( const std::string & file ) 
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
            BarcodeId::Singleton().AssignTag(pair.at(0), std::stoi(pair.at(1)));
        }
        if( in )
            delete in ;
    }
    void BarcodeIdHelper::Print(const std::string & file)
    {
        auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
        if( out == NULL )
            FATAL( " open (barcodeList) file for write failed !!!");
        (*out)<<"0_0_0"<<'\t'<<0<<'\n';
        BarcodeId::Singleton().Print(*out);
        delete out;
    }
}
}
