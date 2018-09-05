#include "common/test/Test.h"
#include "biocommon/fasta/fasta.h"

#include <sstream>

TEST_MODULE_INIT(FastaTest)

using namespace BGIQD::FASTA;

typedef Fasta<SOAP2ContigHead> FastaTest;

typedef FastaReader<FastaTest> Reader;

    std::string  input = ">3 length 64 cvg_0.0_tip_0\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTAGCCAAGCGTG\n\
AGTCAA\n\
>5 length 64 cvg_0.0_tip_1\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTGTCAGGCCCAGT\n\
AGTCATT\n\
>7 length 64 cvg_0.0_tip_0\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAGATGGGGCCTCAG\n\
AGTCATTC\n\
>9 length 64 cvg_0.0_tip_0\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGAGGGAGGTGGGGGAGAA\n\
AGTCATTCA\n\
>11 length 64 cvg_0.0_tip_0\n\
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCAGCTGGGCATGGTGGCGC\n\
AGTCATTCAA\n\
";



TEST(LoadTest01)
{
    std::istringstream ist(input);
    FastaTest test ;

    int i = 0 ;
    while( Reader::LoadNextFasta(ist,test) )
    {
        i++;
        CHECK(unsigned(i*2+1),test.head.contigId);
        CHECK(69+i , test.seq.Len() );
    }
    CHECK( 5 , i );
}

TEST(LoadTest02)
{
    std::istringstream ist(input);
    std::vector<FastaTest> buffer;
    Reader::LoadAllFasta(ist,buffer);
    CHECK( 5 , buffer.size() );
    for( int i = 0 ; i < 5 ; i++ )
    {
        CHECK( 70 + i , buffer[i].seq.Len() );
    }
}
