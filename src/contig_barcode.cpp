#include "contig_barcode.h"
#include "log.h"
#include "logfilter.h"
#include "file_reader.h"
#include "file_writer.h"
#include "sam_parser.h"
#include "stringtools.h"
#include <cassert>
#include <algorithm>
#include <iostream>

namespace BGIQD {
namespace JOB01 {

    using namespace STRING;
    using namespace FILES;
    using namespace SAM;
    using namespace LOG;

    static BarcodeNum bn;
    static logger log1;

    void initLog(const std::string & module)
    {
        logfilter::singleton().get(module,loglevel::INFO,log1);
    }

    int BarcodeNum::barcode2num( const std::string & str )
    {
        std::string key = trim(str);
        auto itr = data.find( key );
        if( itr == data.end() )
        {
            data[key] = next ++ ;
            return data[key];
        }
        return itr->second;
    }
    void BarcodeNum::save( const std::string & file) const
    {
        auto out = FileWriterFactory::GenerateWriterFromFileName(file);
        for( const auto & i : data )
        {
            (*out)<<i.first<<"\t"<<i.second<<std::endl;
        }
        delete out;
    }
    /************************************************************/

    void loadRefBarcodeUniqueInfo( const std::string & file, refBarcodeUniqueInfo & data )
    {
        timer t( log1, std::string("loadRefBarcodeUniqueInfo"));
        auto in = FileReaderFactory::GenerateReaderFromFileName(file);
        while( ! in->eof() )
        {
            std::string line;
            std::getline(*in,line);
            if( in->eof() )
                break;
            auto d1 = split(line," ");
            assert(d1.size() == 2);
            int pos = std::stoi(d1[0]);
            auto d2 = split(d1[1] , "|");
            if( d2.size() < 2 )
            {
                assert( d2[0] == "0 " );
                continue;
            }
            data[pos].insert("0");
            for( size_t i = 1 ; i < d2.size() ; i++ )
            {
                data[pos].insert(trim(d2[i]));
            }
        }
        delete in ;
    }
    void saveRefBarcodeUniqueInfo(const std::string & file , const refBarcodeUniqueInfo & data)
    {
        timer t( log1, std::string("saveRefBarcodeUniqueInfo"));
        auto out = FileWriterFactory::GenerateWriterFromFileName(file);
        for( const auto & pair : data )
        {
            (*out)<<pair.first<<"\t";
            for(auto itr = pair.second.begin() ; itr != pair.second.end() ; )
            {
                (*out<<*itr);
                itr = std::next(itr);
                if ( itr != pair.second.end() )
                {
                    (*out) <<"|";
                }
            }
            (*out)<<std::endl;
        }
        delete out ;
    }
    /************************************************************/

    void loadRefBarcodeInfo( const std::string & file, refBarcodeInfo & data )
    {
        timer t( log1, std::string("loadRefBarcodeInfo"));
        auto in = FileReaderFactory::GenerateReaderFromFileName(file);
        while( ! in->eof() )
        {
            std::string line;
            std::getline(*in,line);
            if( in->eof() )
                break;
            auto d1 = split(line,"\t");
            assert(d1.size() == 2);
            int pos = std::stoi(d1[0]);
            auto d2 = split(d1[1] , "|");
            if( d2.size() < 2 )
            {
                assert( d2[0] == "0 " );
                continue;
            }
            for( size_t i = 1 ; i < d2.size() ; i++ )
            {
                data[pos].push_back(bn.barcode2num(d2[i]));
            }
        }
        bn.save("barcodeId.txt");
        delete in ;
    }

    /************************************************************/

    void loadRefContigInfo( const std::string & file , refContigInfo & data)
    {
        timer t( log1, std::string("loadRefBarcodeInfo"));
        auto in = FileReaderFactory::GenerateReaderFromFileName(file);
        while( ! in->eof() )
        {
            std::string line;
            std::getline(*in,line);
            if( in->eof() )
                break;
            LineParser p(line);
            if(p.IsHead())
                continue;
            auto d0 = p.ParseAsMatchData();
            for( size_t i = 0 ; i< d0.detail.infos.size() ; i++ )
            {
                auto info = d0.detail.infos[i];
                if( info.type != M )
                {
                    continue;
                }
                assert( info.end_position_on_ref - info.start_position_on_ref
                        == info.end_position_on_read - info.start_position_on_read);
                int length = info.end_position_on_ref - info.start_position_on_ref +1 - 63;
                int read = std::stoi(d0.read_name);
                for( int j = 0; j< length; j++ )
                {
                    data[ info.start_position_on_ref + j ].emplace_back( info.start_position_on_read+j , read );
                }
            }
        }
        delete in ;
    }

    /************************************************************/
    void generateConrigBarcodeInfo( const refBarcodeInfo & i_b ,
                                const refContigInfo & i_c,
                                contigBarcodeInfo & data )
    {
        timer t( log1, std::string("generateConrigBarcodeInfo"));
        for( const auto & bi : i_b )
        {
            const auto itr = i_c.find( bi.first ) ;
            if( itr == i_c.end() )
                continue;
            for( const auto & ci : itr->second )
            {
                int pos , contig;
                std::tie( pos , contig )=  ci ;
                data[contig].emplace_back(pos,bi.second);
            }
        }
        for( auto & i : data )
        {
            std::sort(i.second.begin() , i.second.end());
        }
    }

    void printContigBarcodeInfo( const contigBarcodeInfo & data, const std::string & str )
    {
        auto ost = FileWriterFactory::GenerateWriterFromFileName(str);
        timer t( log1, std::string("printContigBarcodeInfo"));
        for( const auto i : data )
        {
            (*ost)<<i.first<<"\t";
            for( const auto & ii : i.second)
            {
                (*ost)<<std::get<0>(ii)<<":";
                int len = std::get<1>(ii).size();
                int index = 1 ;
                for( const auto & iii : std::get<1>(ii) )
                {
                    (*ost)<<iii;
                    if(index < len )
                    {
                        (*ost)<<"|" ;
                        index ++ ;
                    }
                }
                (*ost)<<"\t";
            }
            (*ost)<<std::endl;
        }
        delete ost;
    }
    void loadContigBarcodeInfo( const std::string & file , contigBarcodeInfo & data)
    {
        timer t( log1, std::string("loadContigBarcodeInfo"));
        auto in = FileReaderFactory::GenerateReaderFromFileName(file);
        while( ! in->eof() )
        {
            std::string line;
            std::getline(*in,line);
            if( in->eof() )
                break;
            auto d1 = split(line,"\t");
            assert(d1.size() > 1);
            int contigId = std::stoi(d1[0]);
            for( size_t i = 1 ; i < d1.size() ; i++ )
            {
                auto d2=split(d1[i],":");
                assert(d2.size() == 2);
                int pos = std::stoi(d2[0]);
                auto d3 = split( d2[1] , "|");
                barcodeList list;
                for( const auto &j : d3)
                {
                    list.push_back(std::stoi(j));
                }
                data[contigId].emplace_back(pos,list);
            }
        }
        delete in;
    }
    
    void generateBinBarcodeInfo(const contigBarcodeInfo & data  , int binSize, binBarcodeInfo & d)
    {
        timer t( log1, std::string("generateBinBarcodeInfo"));
        for(const auto & contig : data)
        {
            int contigId = contig.first;
            for( const auto & posData : contig.second )
            {
                int posId= std::get<0>(posData);
                auto v0 = std::get<1>(posData);
                auto  &m = d[contigId] ;

                auto & v= m[posId/binSize];

                for( auto vv : v0)
                {
                    auto itr = v.find(vv) ;
                    if( itr != v.end() )
                    {
                        itr->second ++;
                    }
                    else
                        v[vv] = 1;
                }
                //v.insert(v.end() , v0.begin() , v0.end());
            }
        }
    }

    void saveBinBarcodeInfo(const std::string & file ,const  binBarcodeInfo &data)
    {
        timer t( log1, std::string("saveRefBarcodeUniqueInfo"));
        auto out = FileWriterFactory::GenerateWriterFromFileName(file);
        for( const auto & pair : data )
        {
            for( const auto & pp : pair.second )
            {
                (*out)<<pair.first<<":"<<pp.first;
                for( auto i : pp.second)
                {
                    (*out)<<"\t"<<i.first<<":"<<i.second;
                }
                (*out)<<std::endl;
            }
        }
    }

}//JOB01
}//BGIQD
