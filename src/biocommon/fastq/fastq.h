#ifndef __BIOCOMMON_FASTQ_FASTQ_H__
#define __BIOCOMMON_FASTQ_FASTQ_H__

#include <string>
#include <cassert>
#include <vector>
#include <iostream>

#include "biocommon/seq/seq.h"
#include "common/flags/flags.h"

namespace BGIQD {
    namespace FASTQ {

        struct NormalHead
        {
            std::string head;

            void Init( const std::string & line )
            {
                head = line ;
            }
            std::string Head() const { return head ; }

            void Reset() { head.clear() ;} 
        };

        struct stLFRHeader
        {

            enum ReadType
            {
                Unknow = 0 ,
                readName_barcodeStr =  1 ,
                readName_barcodeStr_index =  2 ,
                readName_barcodeStr_index_barcodeNum = 3 ,
            } type ;


            int readIndex ; //  1/2/3

            int barcode_num;

            std::string barcode_str ;

            std::string readName;

            void Init( const std::string & line ) ;

            std::string Head() const ;

            void Reset()
            {
                type = Unknow ;
                readName = barcode_str = "" ;
                readIndex = barcode_num = 0 ;
            }
        };

        template<class T>
            struct Fastq
            {
                typedef T Header;

                FLAGS_INT ;

                ADD_A_FLAG(1,UnSet);
                ADD_A_FLAG(2,Set_head);
                ADD_A_FLAG(3,Set_seq);
                ADD_A_FLAG(4,Set_3);
                ADD_A_FLAG(5,Set_quality);

                void Reset() 
                {
                    head.Reset();
                    seq.Reset();
                    quality.Reset();
                    flags = 0 ;
                    Set_UnSet();
                }

                Header head;

                BGIQD::SEQ::seq seq;

                BGIQD::SEQ::seq quality;

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

                void Add3(const std::string & )
                {
                    if( Is_UnSet()
                            || ! Is_Set_head() 
                            || ! Is_Set_seq() )
                    {
                        assert(0);
                    }
                    Set_Set_3();
                }

                void AddQuality( const std::string & line )
                {

                    if( Is_UnSet() || ! Is_Set_3() )
                    {
                        assert(0);
                    }
                    Set_Set_quality();
                    quality.AddPartSeq(line);
                }

                bool Is_Setted() const { 
                    return  (!Is_UnSet() )
                        && Is_Set_head() 
                        && Is_Set_seq()
                        && Is_Set_3() 
                        && Is_Set_quality();
                }

                bool QualityFilled() const 
                {
                    assert(quality.Len() <= seq.Len() );
                    return quality.Len() == seq.Len() ;
                }

            };

        template<class T > 
            struct FastqReader
            {
                typedef T Fastq;

                static bool IsHead(const std::string & line)
                {
                    return ( ! line.empty()) && line[0] == '@' ;
                }
                static bool Is_3(const std::string & line)
                {
                    return ( ! line.empty()) && line[0] == '+' ;
                }

                static void LoadAllFastq( 
                        std::istream & ist 
                        , std::vector<Fastq> & buffer )
                {
                    std::string line ;
                    Fastq fq;
                    fq.Reset();
                    bool seq = true ;
                    while ( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            if(fq.Is_Setted())
                            {
                                buffer.push_back(fq);
                            }
                            fq.Reset();
                            fq.AddHead(line);
                        }
                        else if ( Is_3(line) )
                        {
                            fq.Add3(line);
                            seq = !seq ;
                        }
                        else
                        {
                            if( seq )
                                fq.AddSeq(line);
                            else
                                fq.AddQuality(line);
                        }
                    }
                    if(fq.Is_Setted())
                    {
                        buffer.push_back(fq);
                    }
                    fq.Reset();
                }

                static bool LoadNextFasta(std::istream & ist , Fastq & fq)
                {
                    std::string line ;
                    fq.Reset();
                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            fq.AddHead(line);
                            break ;
                        }
                    }
                    if( ist.eof() || ! fq.Is_Set_head() )
                        return false ;

                    bool seq = true ;
                    bool next_head = false;
                    while( ! std::getline(ist,line).eof() )
                    {
                        if( IsHead(line) )
                        {
                            next_head = true ;
                            break ;
                        }
                        else if ( Is_3(line) )
                        {
                            fq.Add3(line);
                            seq = !seq ;
                        }
                        else
                        {
                            if( seq )
                                fq.AddSeq(line);
                            else
                            {
                                fq.AddQuality(line);
                                if( fq.QualityFilled() )
                                {
                                    break ;
                                }
                            }
                        }
                    }
                    if ( (! ist.eof()) && next_head )
                    {
                        // Put the head line back into istream
                        ist.rdbuf()->sputbackc('\n');
                        for( auto  i = line.rbegin() ; i!= line.rend() ; i++ )
                        {
                            ist.rdbuf()->sputbackc(*i);
                        }
                    }
                    return  fq.Is_Setted();
                }
            };
    }
}

#endif //__BIOCOMMON_FASTQ_FASTQ_H__
