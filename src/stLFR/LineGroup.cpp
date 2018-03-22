#include "stLFR/LineGroup.h"
#include "common/files/file_reader.h"
#include <sstream>
#include <cassert>

namespace BGIQD{
    namespace stLFR{

        std::string ContigRoad::to_string() const
        {
            std::ostringstream ost;
            if( headin )
                ost<<'[';
            else
                ost<<'(';
            ost<<group.front()<<'\t'<<*group.rbegin();
            if( tailin )
                ost<<']';
            else
                ost<<')';
            ost<<'\t'<<linear_length<<'\t';

            for( auto j : group)
                ost<<j<<'\t';
            return ost.str();
        }

        void ContigRoad::init(const std::string & buff)
        {
            unsigned int head , tail ;
            char headin1, tailin1,dot;
            std::istringstream ist(buff);
            ist>>headin1>>head>>dot>>tail>>tailin1>>linear_length;

            int length = linear_length;
            if( headin1 == '(' )
            {
                length ++ ;
                headin = false;
            }
            else
                headin = true;
            if( tailin1 == ')')
            {
                length ++ ;
                tailin = false ;
            }
            else
                tailin = true;

            int left = length ;
            while( left > 0 )
            {
                unsigned int next ;
                ist>>next;
                left -- ;
                group.push_back(next);
            }

        }

        std::pair<unsigned int, unsigned int>  ContigRoad::getLinearStep(int index)
        {
            assert( index > 0 && index < linear_length );
            int start = index;
            if(  headin )
                start --;
            return std::make_pair(group[start] , group[start+1]);
        }


        void ContigRoads::LoadRoads(const std::string &file)
        {

            auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(file);
            std::string line;
            while(!std::getline(*in,line).eof())
            {
                ContigRoad r;
                r.init(line);
                roads.push_back(r);
            }
            delete in ;
        }
    }
}