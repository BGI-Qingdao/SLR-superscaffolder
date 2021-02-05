#ifndef __COMMON_ARGS_ARGS_PARSER_H__
#define __COMMON_ARGS_ARGS_PARSER_H__

#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <cassert>
#include <iostream>
#include <getopt.h>

/**********************************************************
 * 
 * Brief    :
 *          define some macros to 
 *           1. automaticly construct Usage.
 *           2. parse command-line parameters.
 *           3. detect the validataion of paramters.
 * 
 * Example  :
 *
 *          START_PARSE_ARGS
 *              DEFINE_ARG_REQUIRED(int,   xxx1 ," explain xxx here1 ") 
 *              DEFINE_ARG_OPTIONAL(float, xxx2 ," explain xxx here2 ", "0.0f")
 *              ......
 *          END_PARSE_ARGS
 *
 * Details :
 *
 *          support type of paramters:
 *              * int
 *              * bool
 *              * float
 *              * string
 *
 *          notice: 
 *              * previous value will be overwrited by 
 *                  call this parameters multi-times 
 *                  except for vector<string>
 *              * bool optional parameter will be setted as false.
 *              * required parameter unsetted will halt the program and print Usage.
 *
 * ********************************************************/

namespace BGIQD{
    namespace ARGS{
        // define details of a parameter as args_union.
        struct args_union
        {
            enum type 
            {
                is_bool = 0,
                is_string = 1,
                is_int = 2,
                is_float = 4,
            };

            union data 
            {
                std::string *s;
                std::vector<std::string> *vs;
                bool b;
                int i;
                long l;
                float f;
            };

            bool optional;

            bool setted ;

            std::string explain;

            std::string name ;

            std::string default_value ;

            type t;

            data d;

            static std::string get_type( type t)
            {
                if( t == is_bool )
                    return "[ no arg ]";
                if( t == is_int)
                    return "[ int arg ]";
                if( t == is_float)
                    return "[ float arg ]";
                if( t == is_string)
                    return "[ string arg ]";
                return "";
            }

            std::string args_to_string() 
            {
                if( t == is_bool )
                    return d.b ? "true" : "false";
                if( t == is_int)
                    return std::to_string(d.i);
                if( t == is_float)
                    return std::to_string(d.f);
                if( t == is_string)
                    return to_string();
                return "";
            }

            args_union(type ty
                    ,const  std::string & n
                    ,bool o
                    ,const std::string & df
                    ,const std::string & exp) 
                : optional(o)
                , setted(false)
                , explain(exp)
                , name(n)
                , default_value(df)
                , t(ty) 

            {
                if(! optional )
                {
                    d.s = NULL;
                }
                if ( t == is_float )
                    d.f = 0.0f ;
                if ( t == is_string )
                    d.s = new std::string("");
            }

            void set_value( const  char * value , bool df )
            {
                if( t == BGIQD::ARGS::args_union::is_bool  )
                {
                    d.b = ! df ;
                }
                if( t == BGIQD::ARGS::args_union::is_string)
                {
                    d.s = new std::string(value);
                }
                if( t == BGIQD::ARGS::args_union::is_int)
                {
                    d.i = std::stoi(std::string(value));
                }
                if( t == BGIQD::ARGS::args_union::is_float)
                {
                    d.f = std::stod(std::string(value));
                }
            }

            std::string to_string() const 
            {
                assert( t == is_string && d.s );
                return *d.s;
            }

            bool to_bool() const { assert(t == is_bool); return d.b ; }

            int to_int() const { assert( t == is_int) ; return d.i ; }


            float to_float() const { assert( t== is_float) ; return d.f ; }

            ~args_union()
            {
                if( t == is_string && d.s )
                    delete d.s ;
            }
        };

        template<class T>
            struct args_traits
            {
                args_union::type type() ;
            };

        //below partial specialization for detail types : 
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
            struct args_traits<float>
            {
                args_union::type type() { return args_union::type::is_float; }
            };

        template<>
            struct args_traits< bool >
            {
                args_union::type type() { return args_union::type::is_bool; }
            };


        // storage for parameters.
        // notice : if those variable names also defined in application file, error may happen.
        static std::map<int,args_union*>  infos;

        static int arg_index = 0 ;

        const int arg_max = 100 ;

        static struct option long_options[arg_max];

    }//ARGS
}//BGIQD

// MARCO for automaticly constructing usage:
#define __PRINT_USAGE \
    std::cerr<<"Usage : "<<argv[0]<< " args "<<std::endl;\
    for( const auto &i : BGIQD::ARGS::infos )\
{\
    std::cerr<<"\t\t"<<"--"<<i.second->name;\
    if( i.second->optional){\
        std::cerr<<"\t"<<"[optional]";}\
    else{std::cerr<<"\t"<<"[required]";}\
    std::cerr<<"\t"<<BGIQD::ARGS::args_union::get_type(i.second->t);\
    std::cerr<<"\t"<<i.second->explain;\
    if( i.second->optional ){\
        std::cerr<<"\t [ default= "<<i.second->default_value<<" ]";\
    }\
    std::cerr<<std::endl;\
}

// MARCO for constructing long optional parameters by default value.
#define __CONSTRUCT_LONG_OPTIONS\
    int max = 0 ;\
    for( const auto &i : BGIQD::ARGS::infos )\
{\
    int id = i.first;\
    auto & item = *(i.second);\
    BGIQD::ARGS::long_options[id].name = item.name.c_str();\
    if(item.t == BGIQD::ARGS::args_union::is_bool )\
    {\
        BGIQD::ARGS::long_options[id].has_arg = 0;\
    }else{\
        BGIQD::ARGS::long_options[id].has_arg = 1;\
    }\
    BGIQD::ARGS::long_options[id].flag = 0 ;\
    BGIQD::ARGS::long_options[id].val = id ;\
    max = id ;\
}\
    for(int i = max + 1 ; i < BGIQD::ARGS::arg_max ; i++ )\
{\
    BGIQD::ARGS::long_options[i].name = 0;\
    BGIQD::ARGS::long_options[i].has_arg = 0;\
    BGIQD::ARGS::long_options[i].flag = 0;\
    BGIQD::ARGS::long_options[i].val = 0;\
}

// MARCO for parsing parameters one by one:
#define __PARSE_ARGS \
    int __curr_flag , __out = 0;\
    while( ( __curr_flag = getopt_long_only( argc , argv,"",BGIQD::ARGS::long_options, &__out  ) ) != EOF )\
{\
    auto itr = BGIQD::ARGS::infos.find( __curr_flag ) ;\
    if ( itr == BGIQD::ARGS::infos.end() )\
    {    continue; }\
    itr->second->setted = true ;\
    itr->second->set_value(optarg, false);\
}\
// MRCRO to print usage
#define __PRINT_ARGS \
    std::cerr<<argv[0];\
    for( const auto &i : BGIQD::ARGS::infos ){\
        if( i.second->setted) {\
            std::cerr<<" --"<<i.second->name<<" "<<i.second->args_to_string();\
        }\
    }\
    std::cerr<<std::endl;

// MRCRO to detect whether all required paremeters are all setted.
#define __CHECK_ARGS \
    bool pass = true ;\
    for( const auto &i : BGIQD::ARGS::infos )\
{\
    if( ! i.second->optional &&  ! i.second->setted ){\
        std::cerr<<"ERROR:  unset nacessary args -- "<<i.second->name<<std::endl;\
        pass =false ;\
    }\
    if( i.second->optional && (! i.second->setted) ){\
        i.second->set_value(i.second->default_value.c_str(),true); \
    }\
}\
if( ! pass ){\
    __PRINT_USAGE\
    return 0;\
} else {\
    __PRINT_ARGS \
}

// construct detail of parameters
#define __DEFINE_ARG_DETAIL( typen , name , optional , d ,exp ) \
    BGIQD::ARGS::args_union name(BGIQD::ARGS::args_traits<typen>().type(),#name, optional,d,exp);\
    BGIQD::ARGS::infos[BGIQD::ARGS::arg_index]=&name;\
    BGIQD::ARGS::arg_index ++ \
// reset storage
#define START_PARSE_ARGS \
    BGIQD::ARGS::infos.clear();\
    BGIQD::ARGS::arg_index = 0 ;\

// wrap for required paramters
#define DEFINE_ARG_REQUIRED( typen , name , exp ) \
    __DEFINE_ARG_DETAIL( typen , name , false,"", exp)

// wrap for optional paramters
#define DEFINE_ARG_OPTIONAL(typen , name ,  exp , df)\
    __DEFINE_ARG_DETAIL( typen , name , true , df , exp);\

// check help command
#define __CHECK_HELP \
    if( argc < 2\
            || std::string(argv[1]) == "-h" \
            || std::string(argv[1]) == "--help"\
            || std::string(argv[1]) == "help"\
            ){\
        __PRINT_USAGE\
        return 0 ;\
    }
// hook all actions here
#define END_PARSE_ARGS \
    __CHECK_HELP\
    __CONSTRUCT_LONG_OPTIONS\
    __PARSE_ARGS\
    __CHECK_ARGS

#endif //__COMMON_ARGS_ARGS_PARSER_H__
