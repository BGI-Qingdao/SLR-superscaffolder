#ifndef __BIOCOMMON_SEQ_SEQ_H__
#define __BIOCOMMON_SEQ_SEQ_H__

#include <string>

namespace BGIQD {
    namespace SEQ {

        struct seq
        {
            std::string atcgs;

            std::string Seq(int weight = -1 ) const ;

            void ToUpper() ;

            void ToLower() ;

            int Len() const ;

            int ValidLen() const ;

            int NLen() const ;

            int GCLen() const ;

            void AddPartSeq( const std::string & line ) ;

        };

    }
}

#endif //__BIOCOMMON_SEQ_SEQ_H__
