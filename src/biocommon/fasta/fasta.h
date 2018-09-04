#ifndef __BIOCOMMON_FASTA_FASTA_H__
#define __BIOCOMMON_FASTA_FASTA_H__

#include <string>
#include <vector>
#include "biocommon/seq/seq.h"

namespace BGIQD {

    namespace FASTA {

        struct IFastaHeader
        {
            virtual void Init( const std::string & line ) = 0;

            virtual std::string Head() const = 0;

        };

        struct IFasta
        {
            
        };

        struct Fasta
        {
            static bool IsHead(const std::string & line) ;

            std::string head;
            std::string seq;

            void SetHead(const std::string & head );
            void AddSeq(const std::string & seq);

            std::string Head() const ;

            std::string Seq(int weight = -1 ) const ;

            int Length() const ;

            //void Seq2Upper() ;
            //void Seq2Lower() ;
        };


        class 

        template<class T = Fasta> 
            struct FastaReader
            {
                std::vector<Fasta> buffer;
            };


    }
}

#endif //__BIOCOMMON_FASTA_FASTA_H__
