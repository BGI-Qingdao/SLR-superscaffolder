#include "stLFR/CBB.h"
#include "sstream"
#include "common/string/stringtools.h"
#include "common/files/file_writer.h"
#include "common/files/file_reader.h"

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
            for( const auto & i : data )
            {
                if( i.empty()) 
                    continue;
                (*out)<<i.ToString()<<std::endl;
            }
            delete out;
        }
    }
}
