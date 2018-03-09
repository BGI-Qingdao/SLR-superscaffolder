#include "soap2/contigGraph.h"
#include "soap2/loadGraph.h"
#include <iostream>
#include "common/args/argsparser.h"
#include "common/files/file_reader.h"
#include "common/string/stringtools.h"
#include "common/log/log.h"
#include "common/log/logfilter.h"
#include "common/multithread/MultiThread.h"
#include "common/freq/freq.h"
BGIQD::LOG::logger lger;

struct A
{
    unsigned int  contig;
    int start ;
    int end;
};
int main(int argc , char **argv)
{
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , lger);
    BGIQD::LOG::timer t(lger,"SuperContig");
    START_PARSE_ARGS
    DEFINE_ARG_DETAIL(std::string , pos, 'p',false,"postion order of contig");
    END_PARSE_ARGS

    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(pos.to_string());
    std::string line;

    std::vector<A> ref_keys;
    std::map<unsigned int,int> ref_keys_map;
    int index = 0 ;
    while(!std::getline(*in,line).eof())
    {
        auto items = BGIQD::STRING::split(line,"\t");
        assert(items.size()>0);
        A a;
        a.contig = std::stoi(items[0]);
        a.start = std::stoi(items[2]);
        a.end= std::stoi(items[3]);
        ref_keys.push_back(a);
        ref_keys_map[std::stoi(items[0])]=index++;
    }
    delete in;

    BGIQD::FREQ::Freq<int>  correct;
    BGIQD::FREQ::Freq<int>  wrong;
    BGIQD::FREQ::Freq<int>  len;
    BGIQD::FREQ::Freq<int>  del;

    int total = 0;
    while(!std::getline(std::cin,line).eof())
    {
        unsigned int head , tail ;
        char headin, tailin,dot;
        int length;
        std::istringstream ist(line);
        ist>>headin>>head>>dot>>tail>>tailin>>length;
        if( length < 2 )
            continue;

        if( headin == '(' )
            length ++ ;
        if( tailin == ')')
            length ++ ;

        total ++ ;
        std::vector<unsigned int> line;
        int left = length ;
        while( left > 0 )
        {
            unsigned int next ;
            ist>>next;
            left -- ;
            if(left == length-1 && headin == '(')
                continue;
            if( left == 0 && tailin == ')')
                continue;
            line.push_back(next);
        }

        auto findR = [&ref_keys, &ref_keys_map]( unsigned int from , unsigned int to )
        {
            int root = ref_keys_map[from];
            //forward
            for( int j = root+1 ; j<(root+4)&& j<(int)ref_keys.size() ; j++)
            {
                if( ref_keys[j].contig == to )
                    return j-root;
            }
            //backward
            for( int j = root-1 ; j>int(root-4)&& j>=0 ; j--)
            {
                if( ref_keys[j].contig == to )
                    return j-root;
            }
            return 0;
        };

        int order = 0;
        bool correct1 = true;
        for( int i =1 ; i< (int)line.size() ; i++ )
        {
            int j = findR(line[i-1],line[i]);
            if( j == 0 )
            {
                wrong.Touch(line.size());
                correct1 = false;
                break;
            }

            if( order == 0)
                order = j;
            else
            {
                // error
                if( order *j <= 0 )
                {
                    wrong.Touch(line.size());
                    correct1 = false;
                    break;
                }
                order += j ;
            }
        }
        if( !correct1)
            continue;
        bool go = order > 0 ? true : false;
        if ( order < 0 ) order = -order;
        if( order == (int)(line.size()) -1 )
        {
            correct.Touch(line.size());
            if( go )
                len.Touch(std::abs( ref_keys[ref_keys_map[line[0]]].start- ref_keys[ref_keys_map[*line.rbegin()]].end ));
            else
                len.Touch(std::abs( ref_keys[ref_keys_map[line[0]]].end- ref_keys[ref_keys_map[*line.rbegin()]].start));
        }
        else
        {
            del.Touch(order - line.size() +1 );
        }
    }
    std::cout<<"correct freq\n"<<correct.ToString()<<std::endl;
    std::cout<<"wrong freq\n"<<wrong.ToString()<<std::endl;
    std::cout<<"del freq\n"<<del.ToString()<<std::endl;
    std::cout<<"len freq\n"<<len.ToString()<<std::endl;
}
