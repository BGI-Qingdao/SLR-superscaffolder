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
    DEFINE_ARG_REQUIRED(std::string , pos, "postion order of contig");
    END_PARSE_ARGS

    auto in = BGIQD::FILES::FileReaderFactory::GenerateReaderFromFileName(pos.to_string());
    std::string line_1;

    std::vector<A> ref_keys;
    std::map<unsigned int,std::vector<int> > ref_keys_map;
    int index = 0 ;
    while(!std::getline(*in,line_1).eof())
    {
        auto items = BGIQD::STRING::split(line_1,"\t");
        assert(items.size()>0);
        A a;
        a.contig = std::stoi(items[0]);
        a.start = std::stoi(items[2]);
        a.end= std::stoi(items[3]);
        ref_keys.push_back(a);
        ref_keys_map[std::stoi(items[0])].push_back(index++);
    }
    delete in;

    BGIQD::FREQ::Freq<int>  correct;
    BGIQD::FREQ::Freq<int>  wrong;
    BGIQD::FREQ::Freq<int>  len;
    BGIQD::FREQ::Freq<int>  del;
    BGIQD::FREQ::Freq<unsigned int>  seeds;

    int total = 0;
    while(!std::getline(std::cin,line_1).eof())
    {
        unsigned int head , tail ;
        char headin, tailin,dot,dir;
        int length;
        std::istringstream ist(line_1);
        ist>>dir>>headin>>head>>dot>>tail>>tailin>>length;
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
        /*
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
        */
        for( int i =0 ; i< (int)line.size() ; i++ )
        {
            seeds.Touch(line[i]);
        }

        int status = 0 ;
        int del_total = 0 ;
        for( const auto root : ref_keys_map[line[0]])
        {
            //downstream
            int order = 1 ;
            int start = root + order ;
            int i = 1 ;
            del_total = 0 ;
            for( i = 1 ; i< (int)line.size() ; i++ )
            {
                int del = 0;
                while(1)
                {
                    if( ref_keys[ start ].contig == line[i] )
                    {
                        start += order;
                        break;
                    }
                    else
                    {
                        del ++ ;
                        start += order;
                    }
                    if( del > 3 )
                        break;
                }
                if( del > 3 ) 
                    break;
                else
                    del_total += del;
            }
            if( i == (int)line.size() )
            {
                status = 1 ; //correct
                break;
            }
            //upstream
            order = -1 ;
            start = root + order ;
            del_total = 0 ;
            i = 1 ;
            for( i = 1 ; i< (int)line.size() ; i++ )
            {
                int del = 0;
                while(1)
                {
                    if( ref_keys[ start ].contig == line[i] )
                    {
                        start += order;
                        break;
                    }
                    else
                    {
                        del ++ ;
                        start += order;
                    }
                    if( del > 5 )
                        break;
                }
                if( del > 5 ) 
                    break;
                else
                    del_total +=del;
            }
            if( i == (int)line.size() )
            {
                status = 1 ; //correct
                break;
            }
        }
        //int del 
        /*
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
        */
        if( status == 1 && del_total == 0 )
            correct.Touch(line.size());
        else if ( status == 1 && del_total > 0 )
        {
            del.Touch(del_total);
            std::cout<<"D: "<<line_1<<std::endl;
        }
        else
        {
            wrong.Touch(line.size());
            std::cout<<"E: "<<line_1<<std::endl;
        }
    }
    std::cout<<"correct freq\n"<<correct.ToString()<<std::endl;
    std::cout<<"wrong freq\n"<<wrong.ToString()<<std::endl;
    std::cout<<"del freq\n"<<del.ToString()<<std::endl;
    //std::cout<<"len freq\n"<<len.ToString()<<std::endl;
    for( const auto i: seeds.data)
    {
        if( i.second > 1 )
        {
            std::cout<<"err 2 "<<i.first<<"\t"<<i.second<<std::endl;
        }
    }
}
