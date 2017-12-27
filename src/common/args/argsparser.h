#ifndef __COMMON_ARGS_ARGS_PARSER_H__
#define __COMMON_ARGS_ARGS_PARSER_H__

#include <string>
#include <tuple>
#include <vector>
#include <stdarg.h>
#include <unistd.h>
#include <map>

namespace BGIQD{
namespace ARGS{

struct args_union
{
    enum type
    {
        is_bool = 0,
        is_string = 1,
        is_int = 2,
        is_long = 3,
        is_float = 4,
        is_vector_string = 5,
    };

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

    args_union() : t(is_bool) {}

    std::string to_string() const 
    {
        return *d.s;
    }
    std::vector<std::string> to_vector_string() const 
    {
        return *d.vs;
    }

    bool to_bool() const { return d.b ; }
    int to_int() const { return d.i ; }
    long to_long() const { return d.l ; }
    float to_float() const { return d.f ; }
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

#define START_PARSE_ARGS  

#define DEFINE_ARG( typen , name , flag ) \
    args_union name;\
    name.t= args_traits<typen>().type();\
    infos[flag]=&name;


#define END_PARSE_ARGS \
    std::string format;\
    for( const auto &i : infos )\
    {\
        format += i.first;\
        if(i.second->t != args_union::is_bool )\
        {\
            format += ":";\
        }\
    }\
    int curr_flag ;\
    while( ( curr_flag = getopt( argc , argv, format.c_str() ) ) != EOF )\
    {\
        auto itr = infos.find( curr_flag ) ;\
        if ( itr == infos.end() )\
        {    continue; }\
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


#endif //__COMMON_ARGS_ARGS_PARSER_H__
