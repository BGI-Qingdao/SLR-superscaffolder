#include "utils/misc/contigIndex.h"

#include <sstream>
#include <cassert>

namespace BGIQD {
    namespace SOAP2 {

        int ContigIndex::K  ;

        std::string ContigIndex::ToString() const 
        {
            std::ostringstream ist;
            ist<<contig<<'\t'<<length<<'\t'<<reverse_add;
            return ist.str();
        }
        void ContigIndex::InitFromString(const std::string &line )
        {
            std::istringstream ist(line);
            ist>>contig>>length;
            if( ! ist.eof() )
                ist>>reverse_add;
            else
                reverse_add = 1 ;
        }


        void ContigIndexMap::LoadContigIndexs( std::istream &ist )
        {
            std::string line ;
            ContigIndex tmp ;

            while( ! std::getline(ist,line).eof())
            {
                if(! std::isdigit( line[0]) )
                    continue ;
                tmp.InitFromString(line);
                data[tmp.contig] = tmp ;
            }
        }

        void ContigIndexMap::BuildReverseCompleteContigs()
        {
            std::map<unsigned int , unsigned int > c2;
            for( const auto &pair : data)
            {
                const auto & c1 = pair.second ;
                if( c1.reverse_add == 1 )
                {
                    assert( data.find(c1.contig+1 ) == data.end());
                    c2[c1.contig+1] = c1.contig ;
                }
            }
            for( const auto & pair : c2 )
            {
                try
                {
                    data[pair.first] = data.at(pair.second) ;
                }
                catch(...)
                {
                    assert(0);
                }
            }
        }

        const ContigIndex & ContigIndexMap::GetContigIndex( unsigned int id ) const 
        {
            try
            {
                return data.at(id) ;
            }
            catch(...)
            {
                assert(0);
            }
            static ContigIndex tmp;
            // never come here .
            return tmp;
        }
    }
}
