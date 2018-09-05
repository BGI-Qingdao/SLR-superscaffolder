#ifndef __BIOCOMMON_FASTA_FASTA_H__
#define __BIOCOMMON_FASTA_FASTA_H__

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#include "biocommon/seq/seq.h"
#include "common/flags/flags.h"

namespace BGIQD {

    namespace FASTA {

        struct IFastaHeader
        {
            virtual void Init( const std::string & line ) = 0;

            virtual std::string Head() const = 0;

            virtual void Reset() = 0 ;
        };

        struct NormalHead : public IFastaHeader
        {
            std::string head;

            virtual void Init( const std::string & line ) final
            {
                head = line ;
            }
            virtual std::string Head() const { return head ; }

            virtual void Reset() final { head.clear() ;} 
        };

        struct SOAP2ContigHead : public IFastaHeader
        {
            unsigned int contigId ;

            int len ;

            float cov;

            int is_tip ;

            virtual void Init( const std::string & line ) final;

            virtual std::string Head() const final ;

            virtual void Reset() final { contigId = 0 ; is_tip = 0 ; cov = 0 ; len = 0;  } ;
        };

        struct SOAP2ScaffHead : public IFastaHeader 
        {
            std::string head;

            virtual void Init( const std::string & line ) final
            {
                head = line ;
            }
            virtual std::string Head() const { return head ; }

            virtual void Reset() final { head.clear() ;} 

        };

        template<class T>
            struct Fasta
            {
                typedef T Header;

                FLAGS_INT ;

                ADD_A_FLAG(1,UnSet);
                ADD_A_FLAG(2,Set_head);
                ADD_A_FLAG(3,Set_seq);

                void Reset() 
                {
                    head.Reset();
                    seq.Reset();
                    flags = 0 ;
                    Set_UnSet();
                }

                Header head;
                BGIQD::SEQ::seq seq;
                void AddHead(const std::string & line)
                {
                    Clean_UnSet();
                    Set_Set_head();
                    head.Init(line);
                }

                void AddSeq(const std::string & line )
                {
                    if( Is_UnSet() || ! Is_Set_head() )
                    {
                        assert(0);
                    }
                    Set_Set_seq();
                    seq.AddPartSeq(line);
                }
                bool Is_Setted() const { 
                    return  (!Is_UnSet() )
                        && Is_Set_head() 
                        && Is_Set_seq() ;
                }

            };

        template<class T > 
            struct FastaReader
            {
                typedef T Fasta;

                static bool IsHead(const std::string & line)
                {
                    return (! line.empty()) && line[0] == '<' ;
                }

                static void LoadAllFasta( std::istream & ist , std::vector<Fasta> & buffer )
                {
                    std::string line ;
                    Fasta fa;
                    fa.Reset();
                    while ( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            if(fa.Is_Setted())
                            {
                                buffer.push_back(fa);
                            }
                            fa.Reset();
                            fa.AddHead(line);
                        }
                        else
                        {
                            fa.AddSeq(line);
                        }
                    }
                    if(fa.Is_Setted())
                    {
                        buffer.push_back(fa);
                    }
                    fa.Reset();
                }

                static bool LoadNextFasta(std::istream & ist , Fasta & fa)
                {
                    std::string line ;

                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            fa.head.Init(line);
                            break ;
                        }
                    }
                    if( ist.eof() )
                        return false ;

                    bool seq_set = false ;

                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            fa.head.Init(line);
                            break ;
                        }
                        else
                        {
                            seq_set = true ;
                            fa.seq.AddSeq(line);
                        }
                    }
                    if (! ist.eof() )
                    {
                        ist.rdbuf()->sputc('\n');
                        ist.rdbuf()->sputn(line.c_str(),line.size());
                    }
                    return  seq_set ;
                }
            };
    }
}

#endif //__BIOCOMMON_FASTA_FASTA_H__
