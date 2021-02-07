#ifndef __COMMON_LOG_LOG_H__
#define __COMMON_LOG_LOG_H__

#include <string>
#include <iostream>
#include <sys/time.h>
#include <sstream>

/**********************************************************
 *
 * @Brief :
 *
 *   * To make uniform logs by easy-to-use interface.
 *   * Automatic benchmark the cpu and wall clock time.
 *
 * *******************************************************/

namespace BGIQD{
    namespace LOG{

        // empty class for template
        struct lstart {};
        struct lend{};

        // stream like logger :
        // logger<<lstart()<<xxx<<xxx<<yyy<<lend();
        class logger
        {
            public:
                void Init(std::string name) 
                {
                    module = name ;
                }
                // clean buffer
                logger & operator << (const lstart & )
                {
                    buffer.str("");
                    return *this;
                }
                // Cache all log
                template< class T >
                    logger & operator << (const T & t)
                    {
                        buffer<<t;
                        return *this;
                    }
                // print all log into stderr
                logger & operator << (const lend & )
                {
                    std::cerr<<logstring(buffer.str())<<std::endl;
                    buffer.str("");
                    return *this;
                }

            private:
                std::ostringstream buffer;
                std::string logstring(const std::string & str);
                std::string module;
        };

        struct timeperoid;

        // A special time point. 
        struct timepoint
        {
            public:
                timeval wall;
                clock_t cpu;

            public:
                static timepoint now() 
                {
                    timepoint ret ;
                    gettimeofday( &ret.wall ,NULL);
                    ret.cpu = clock();
                    return ret;
                }

                timeperoid operator-(const timepoint & prev_point ) const ;
                std::string to_string() const ;
        };

        // A special time period. 
        struct timeperoid 
        {
            public:
                timeval wall;
                clock_t cpu;
            public:
                std::string to_string() const ;
        };

        // timer to log the running time use above structures.
        class timer
        {
            public:
                // get a time tag while constructing
                timer( logger & a_logger , const std::string &job_description)
                    : start(timepoint::now())
                      , l(a_logger)
                      , jobdec(job_description) {
                          (l)<<lstart()<<jobdec<< " start now ... "<<lend();
                      }
                // calculate time period and print it when destructing.
                ~timer()
                {
                    timepoint end = timepoint::now();
                    std::string last = (end-start).to_string();
                    (l)<<lstart()<<jobdec<< " finish. used "<<last<<lend();
                }
            private:
                timepoint start ;
                logger & l;
                std::string jobdec;
        };

    }//LOG
}//BGIQD

#endif //__COMMON_LOG_LOG_H__
