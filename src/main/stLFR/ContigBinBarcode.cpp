#include "stLFR/ContigBinBarcode.h"
#include "sstream"
#include "utils/string/stringtools.h"
#include "utils/files/file_writer.h"
#include "utils/files/file_reader.h"
#include "utils/misc/Error.h"

namespace BGIQD{
    namespace stLFR {

        bool BarcodeOnBin::empty() const 
        {
            return collections.size() == 0 ;
        }

        std::string BarcodeOnBin::ToString() const 
        {
            std::ostringstream ost;
            ost<<contigId<<':'<<binId<<':'<<start<<':'<<end;
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
            start =  std::stoi(d0[2]);
            end =  std::stoi(d0[3]);

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
            ost<<contigId<<':'<<binId<<':'<<start<<':'<<end;
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
            ist>>contigId>>split>>binId>>split>>start>>split>>end;
            int i = 0;
            while( ! ist.eof() )
            {
                BinSimularity sinfo;
                ist>>sinfo.contigId>>split
                    >>sinfo.binId>>split
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

        std::string ContigBarcodeInfo::ToString() const 
        {
            std::ostringstream ost;
            ost<<contig_id<<'\t';
            for( const auto & i : barcodesOnPos)
            {
                bool first = true ;
                ost<<i.first<<":";
                for( unsigned int barcode : i.second )
                {
                    if( first )
                        first = false ;
                    else
                        ost<<'|';
                    ost<<barcode;
                }
                ost<<'\t';
            }
            return ost.str();
        }

        void ContigBarcodeInfo::InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            ist>>contig_id;
            std::string tmp;
            while(!ist.eof())
            {
                ist>>tmp;
                auto item1 = BGIQD::STRING::split(tmp,":");
                int pos = std::stoi(item1[0]);
                auto item2 = BGIQD::STRING::split(item1[1],"|");
                for(int i = 0 ; i < (int)item2.size() ; i++ )
                {
                    barcodesOnPos[pos].insert(std::stoi(item2[i]));
                }
            }
        }

        std::string ContigOnBarcode::ToString() const 
        {
            std::ostringstream ost ;
            ost<<barcode_id;
            for(const auto & pair : contig_data )
            {
                ost<<'\t'<<pair.first<<'\t'<<pair.second ;
            }
            return ost.str();
        }

        void ContigOnBarcode::InitFromString(const std::string & line)
        {
            std::istringstream ist(line);
            ist>>barcode_id;
            while(!ist.eof())
            {
                unsigned int contig;
                int num ;
                ist>>contig>>num ;
                contig_data[contig] = num ;
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
                if( data.at(i).empty() )
                    continue ;
                (*out)<<data.at(i).ToString()<<std::endl;
            }
            delete out;
        }
    }
}
