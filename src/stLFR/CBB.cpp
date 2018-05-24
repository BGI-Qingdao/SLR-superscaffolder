#include "stLFR/CBB.h"
#include "sstream"
#include "common/string/stringtools.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"
#include "common/error/Error.h"

namespace BGIQD{
    namespace stLFR {

        bool BarcodeOnBin::empty() const 
        {
            return collections.size() == 0 ;
        }

        std::string BarcodeOnBin::ToString() const 
        {
            std::ostringstream ost;
            ost<<contigId<<':'<<binId;
            for( auto i : collections)
            {
                ost<<'\t'<<i.first<<':'<<i.second;
            }
            return ost.str();
        }

        void BarcodeOnBin::InitFromString(const std::string & line)
        {
            auto d1 = BGIQD::STRING::split(line,"\t");
            assert(d1.size() >1  );
            auto d0 = BGIQD::STRING::split(d1[0],":");
            contigId = std::stoul(d0[0]);
            binId = std::stoi(d0[1]);
            for( size_t i = 1 ; i < d1.size(); i++)
            {
                auto d2 = BGIQD::STRING::split(d1[i],":");
                int barcodeId = stoi(d2[0]);
                int num = stoi(d2[1]);
                collections.IncreaseElement(barcodeId,num);
            }
        }

        void LoadBarcodeOnBinArray( const std::string & file , BarcodeOnBinArray & data )
        {
            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( " open prefix.barcodeOnBin to read failed !!!" );
            auto add_data = [&] ( const std::string & line)
            {
                BarcodeOnBin b2b ;
                b2b.InitFromString(line);
                data.push_back(b2b);
            };

            BGIQD::FILES::FileReaderFactory::EachLine(*in,add_data);
            delete in ;
        }

        void PrintBarcodeOnBinArray( const std::string & file , const BarcodeOnBinArray & data)
        {
            auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
            if( out == NULL )
                FATAL( " open prefix.barcodeOnBin to write failed !!!" );
            for( const auto & i : data )
            {
                if( i.empty()) 
                    continue;
                (*out)<<i.ToString()<<std::endl;
            }
            delete out;
        }

        std::string BinRelation::ToString() const
        {
            std::ostringstream ost;
            ost<<contigId<<':'<<binId;
            for( auto pair : sims )
            {
                auto & sinfo = pair.second ;
                ost<<'\t'<<sinfo.contigId<<':'<<sinfo.binId<<':'<<sinfo.simularity;
            }
            return ost.str();
        }

        void BinRelation::InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            char split;
            ist>>contigId>>split>>binId;
            int i = 0;
            while( ! ist.eof() )
            {
                BinSimularity sinfo;
                ist>>sinfo.contigId>>split
                    >>sinfo.binIndex>>split
                    >>sinfo.simularity;
                sinfo.binIndex = i ; //Not valid now
                sims[i] = sinfo ;
                i++ ;
            }
        }

        void LoadBinRelationArray(const std::string & file , BinRelationArray & data)
        {
            auto parseline = [&data](const std::string &line)
            {
                BinRelation a_data ;
                a_data.InitFromString(line);
                data.push_back(a_data);
            };

            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( " open prefix.bin_cluster to read failed !!!" );
            BGIQD::FILES::FileReaderFactory::EachLine(*in ,parseline);
            delete in ;
        }

        void PrintBinRelationArray(const std::string & file ,const BinRelationArray & data)
        {
            auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file);
            if( out == NULL )
                FATAL( " open prefix.bin_cluster to write failed !!!" );
            for( size_t i = 0 ; i < data.size() ; i++ )
            {
                (*out)<<data.at(i).ToString()<<std::endl;
            }
            delete out;
        }

        void ContigRelation::InitFromString(const std::string &line )
        {
            std::istringstream ist(line);
            ist>>contigId;

            while(! ist.eof() )
            {
                ContigSimularity sinfo;
                ist>>sinfo.contigId>>sinfo.simularity;
                sims[sinfo.contigId] = sinfo ;
            }
        }

        std::string ContigRelation::ToString() const 
        {
            std::ostringstream ost ;
            ost<<contigId;
            for( const auto pair : sims)
            {
                ost<<'\t'<<pair.second.contigId<<'\t'<<pair.second.simularity;
            }
            return ost.str();
        }

        void LoadContigRelationArray(const std::string & file ,ContigRelationArray & data)
        {
            auto parseline = [&data](const std::string &line)
            {
                ContigRelation a_data ;
                a_data.InitFromString(line);
                data.push_back(a_data);
            };

            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            if( in == NULL )
                FATAL( " open prefix.cluster to read failed !!!" );
            BGIQD::FILES::FileReaderFactory::EachLine(*in ,parseline);
            delete in ;
        }

        void PrintContigRelationArray(const std::string & file ,const ContigRelationArray & data)
        {
            auto out = BGIQD::FILES::FileWriterFactory::GenerateWriterFromFileName(file) ;
            if( out  == NULL )
                FATAL( " open prefix.cluster to write failed !!!" );
            for( size_t i = 0 ; i <data.size() ; i++)
            {
                (*out)<<data.at(i).ToString()<<std::endl;
            }
            delete out;
        }
    }
}
