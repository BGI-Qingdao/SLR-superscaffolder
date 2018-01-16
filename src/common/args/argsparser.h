#ifndef __COMMON_ARGS_ARGS_PARSER_H__
#define __COMMON_ARGS_ARGS_PARSER_H__

#include <string>
#include <tuple>
#include <vector>
#include <stdarg.h>
#include <unistd.h>
#include <map>
#include <cassert>
#include <iostream>

namespace BGIQD{
namespace ARGS{

struct args_union
{
    bool optional;
    bool setted ;
    std::string explain;
    enum type 
    {
        is_bool = 0,
        is_string = 1,
        is_int = 2,
        is_long = 3,
        is_float = 4,
        is_vector_string = 5,
    };

    static std::string get_type( type t)
    {
        if( t == is_bool )
            return " ";
        if( t == is_int)
            return "int";
        if( t == is_long)
            return "long";
        if( t == is_float)
            return "float";
        if( t == is_string)
            return "string";
        if( t == is_vector_string)
            return "vector_string";
        return "";
    }

    std::string args_to_string() 
    {
        if( t == is_bool )
            return d.b ? "true" : "false";
        if( t == is_int)
            return std::to_string(d.i);
        if( t == is_long)
            return std::to_string(d.l);
        if( t == is_float)
            return std::to_string(d.f);
        if( t == is_string)
            return to_string();
        if( t == is_vector_string)
        {
            std::string ret ;
            for(const auto & i : to_vector_string())
            {
                ret += "\t";
                ret += i;
            }
            return ret ;
        }
        return "";
    }

    union data {
        std::string *s;
        std::vector<std::string> *vs;
        bool b;
        int i;
        long l;
        float f;
    };

    type t;
    data d;

    args_union(type t , bool o = true ) : optional(o),setted(false), t(t) {
        if(! optional )
        {
            d.s = NULL;
        }
        if ( t != is_string && t != is_vector_string )
            d.l = 0 ;
        if ( t == is_float )
            d.f = 0.0f ;
        if ( t == is_string )
            d.s = new std::string("");
        if ( t == is_vector_string )
            d.vs = new std::vector<std::string>();
    }

    std::string to_string() const 
    {
        assert( t == is_string && d.s );
        return *d.s;
    }
    std::vector<std::string> to_vector_string() const 
    {
        assert( t == is_vector_string && d.vs );
        return *d.vs;
    }

    bool to_bool() const { assert(t == is_bool); return d.b ; }
    int to_int() const { assert( t == is_int) ; return d.i ; }
    long to_long() const { assert(t == is_long) ;return d.l ; }
    float to_float() const { assert( t== is_float) ; return d.f ; }

    template<class T> 
    T value(T &)
    {
        assert(0);
    }

    template<class T=int>
    int value(int &a)
    {
        a= to_int();
        return a;
    }

    template<class T=long>
    long value(long &a)
    {
        a= to_long();
        return a;
    }
    template<class T=float>
    float value(float &a)
    {
        a= to_float();
        return a;
    }

    template<class T=bool>
    bool value(bool &a)
    {
        a= to_bool();
        return a;
    }

    template<class T=std::string>
    std::string value(std::string &a)
    {
        a= to_string();
        return a;
    }
    template<class T=std::vector<std::string>>
    std::vector<std::string> value(std::vector<std::string> &a)
    {
        a= to_vector_string();
        return a;
    }

    ~args_union()
    {
        if( t == is_string && d.s )
            delete d.s ;
        if ( t== is_vector_string && d.vs )
            delete d.vs ;
    }
};

static std::map<int,args_union*>  infos;

template<class T>
struct args_traits
{
    args_union::type type() ;
};

template<>
struct args_traits<int>
{
    args_union::type type() { return args_union::type::is_int ; }
};

template<>
struct args_traits<std::string>
{
    args_union::type type() { return args_union::type::is_string; }
};

template<>
struct args_traits<long>
{
    args_union::type type() { return args_union::type::is_int ; }
};


template<>
struct args_traits<float>
{
    args_union::type type() { return args_union::type::is_float; }
};

template<>
struct args_traits<std::vector<std::string> >
{
    args_union::type type() { return args_union::type::is_vector_string; }
};

template<>
struct args_traits< bool >
{
    args_union::type type() { return args_union::type::is_bool; }
};


}//ARGS
}//BGIQD

#define __PRINT_USAGE \
    std::cerr<<"Usage : "<<argv[0]<< " args "<<std::endl;\
    for( const auto &i : infos )\
    {\
        std::cerr<<"\t\t"<<"-"<<(char)i.first;\
        if( i.second->optional){\
            std::cerr<<"\t"<<"[optional]";}\
        std::cerr<<"\t"<<BGIQD::ARGS::args_union::get_type(i.second->t)\
                 <<"\t"<<i.second->explain<<std::endl;\
    }
#define __CONSTRUCT_FORMAT \
    std::string format;\
    for( const auto &i : infos )\
    {\
        format += i.first;\
        if(i.second->t != args_union::is_bool )\
        {\
            format += ":";\
        }\
    }\

#define __PARSE_ARGS \
    int curr_flag ;\
    while( ( curr_flag = getopt( argc , argv, format.c_str() ) ) != EOF )\
    {\
        auto itr = infos.find( curr_flag ) ;\
        if ( itr == infos.end() )\
        {    continue; }\
        itr->second->setted = true ;\
        if( itr->second->t == args_union::is_bool )\
        {\
            itr->second->d.b = true;\
        }\
        if( itr->second->t == args_union::is_string)\
        {\
            itr->second->d.s = new std::string(optarg);\
        }\
        if( itr->second->t == args_union::is_int)\
        {\
            itr->second->d.i = std::stoi(std::string(optarg));\
        }\
        if( itr->second->t == args_union::is_long)\
        {\
            itr->second->d.l = std::stol(std::string(optarg));\
        }\
        if( itr->second->t == args_union::is_float)\
        {\
            itr->second->d.f = std::stod(std::string(optarg));\
        }\
        if( itr->second->t == args_union::is_vector_string)\
        {\
            if( itr->second->d.vs == NULL )\
            {\
                itr->second->d.vs = new std::vector<std::string>() ;\
            }\
            (*(itr->second->d.vs)).push_back(std::string(optarg));\
        }\
    }\

#define __PRINT_ARGS \
    std::cerr<<argv[0];\
    for( const auto &i : infos ){\
        if( i.second->setted) {\
        std::cerr<<"\t-"<<(char)i.first<<"\t"<<i.second->args_to_string();\
        }\
    }\
    std::cerr<<std::endl;

#define __CHECK_ARGS \
    bool pass = true ;\
    for( const auto &i : infos )\
    {\
        if( ! i.second->optional &&  ! i.second->setted ){\
            std::cerr<<"ERROR:  unset nacessary args - "<<(char)i.first<<std::endl;\
            pass =false ;\
        }\
    }\
    if( ! pass ){\
        __PRINT_USAGE\
        return 0;\
    } else {\
        __PRINT_ARGS \
    }

#define __DEFINE_ARG_DETAIL( typen , name , flag , optional ) \
    args_union name(args_traits<typen>().type(), optional);\
    infos[flag]=&name;

#define START_PARSE_ARGS  

#define DEFINE_ARG_REQUIRED( typen , name , flag ) \
    __DEFINE_ARG_DETAIL( typen , name , flag , false)

#define DEFINE_ARG( typen , name , flag ) \
    __DEFINE_ARG_DETAIL( typen , name , flag , true )

#define DEFINE_ARG_DETAIL(typen , name , flag ,o,  e)\
    __DEFINE_ARG_DETAIL( typen , name , flag , o);\
    name.explain = e ;

#define END_PARSE_ARGS \
    __CONSTRUCT_FORMAT\
    __PARSE_ARGS\
    __CHECK_ARGS

#endif //__COMMON_ARGS_ARGS_PARSER_H__
