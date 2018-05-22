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
    bool is_reverse_complete;
};

int main(int argc , char **argv)
{
    BGIQD::LOG::logfilter::singleton().get("SuperContig",BGIQD::LOG::loglevel::INFO , lger);
    BGIQD::LOG::timer t(lger,"SuperContig");
    START_PARSE_ARGS
    DEFINE_ARG_REQUIRED(std::string , pos, "postion order of contig");
    DEFINE_ARG_OPTIONAL(bool , detail, "print error detail","0");
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
        a.is_reverse_complete = ( items[4] == "1" );
        ref_keys.push_back(a);
        ref_keys_map[std::stoul(items[0])].push_back(index++);
    }
    delete in;

    BGIQD::FREQ::Freq<int>  correct;
    BGIQD::FREQ::Freq<int>  wrong;
    BGIQD::FREQ::Freq<int>  len;
    BGIQD::FREQ::Freq<int>  del;
    BGIQD::FREQ::Freq<unsigned int>  seeds;

    while(!std::getline(std::cin,line_1).eof())
    {
        auto items = BGIQD::STRING::split(line_1,"\t");

        int length = std::stoi( items[0] ) ;

        if ( length < 2 ) 
            continue;

        std::vector<unsigned int> line;
        std::vector<bool> orders;

        for( int i = 1 ; i <= length ; i++ )
        {
            unsigned int contig = std::stoul(items[i]);
            if( ref_keys_map.find(contig) == ref_keys_map.end() 
                    &&  ref_keys_map.find(contig-1) != ref_keys_map.end())
            {
                line.push_back(contig-1);
                orders.push_back(false);
            }
            else
            {
                line.push_back(contig);
                orders.push_back(true);
            }
        }

        for( int i =0 ; i< (int)line.size() ; i++ )
        {
            seeds.Touch(line[i]);
        }

        int status = 0 ;

        int del_total = 0 ;

        for( const auto root : ref_keys_map[line[0]])
        {
            auto & contig = ref_keys[root] ;
            bool downstream = true ;
            // real contig upstream on ref
            if(  (!contig.is_reverse_complete) && (! orders[0] ) )
            {
                downstream = false ;
            }
            // seed contig upstream on ref
            if ( contig.is_reverse_complete && orders[0] )
            {
                downstream = false ;
            }
            //downstream
            int order = downstream ? 1 : -1 ;
            int start = root + order ;
            int i = 1 ;
            del_total = 0 ;
            bool conflict = false ;
            if ( downstream )
            {
                for( i = 1 ; i< (int)line.size() ; i++ )
                {
                    int del = 0;
                    while(1)
                    {
                        if( ref_keys[ start ].contig == line[i] )
                        {
                            if ( orders[i] && ( ref_keys[start].is_reverse_complete )) 
                            {
                                conflict = true ;
                            }
                            if ( (! orders[i]) && ( !ref_keys[start].is_reverse_complete )) 
                            {
                                conflict = true ;
                            }
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
                    if ( conflict )
                    {
                        break ;
                    }
                    if( del > 3 ) 
                        break;
                    else
                        del_total += del;
                }
                if( ! conflict && i == (int)line.size() )
                {
                    status = 1 ; //correct
                    break;
                }
            }
            else
            {
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
                            if ( orders[i] && ( !ref_keys[start].is_reverse_complete )) 
                            {
                                conflict = true ;
                            }
                            if ( (! orders[i]) && ( ref_keys[start].is_reverse_complete )) 
                            {
                                conflict = true ;
                            }
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
                if( ! conflict && i == (int)line.size() )
                {
                    status = 1 ; //correct
                    break;
                }
            }
        }
        if( status == 1 && del_total == 0 )
            correct.Touch(line.size());
        else if ( status == 1 && del_total > 0 )
        {
            del.Touch(line.size());
            if( detail.to_bool())
                std::cout<<"D: "<<line_1<<std::endl;
        }
        else
        {
            wrong.Touch(line.size());
            if( detail.to_bool())
                std::cout<<"E: "<<line_1<<std::endl;
        }
    }

    std::cout<<"correct freq\n"<<correct.ToString()<<std::endl;
    std::cout<<"wrong freq\n"<<wrong.ToString()<<std::endl;
    std::cout<<"del freq\n"<<del.ToString()<<std::endl;
    //std::cout<<"len freq\n"<<len.ToString()<<std::endl;
    int multi_node = 0 ;
    for( const auto i: seeds.data)
    {
        if( i.second > 1 )
        {
            multi_node ++ ;
            if( detail.to_bool())
                std::cout<<"err 2 "<<i.first<<"\t"<<i.second<<std::endl;
        }
    }

    std::cout<<"multi_node  "<<multi_node<<std::endl; 

    return 0;
}
