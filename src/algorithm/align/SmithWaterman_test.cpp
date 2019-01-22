#include "algorithm/align/SmithWaterman.h"
#include "common/test/Test.h"

TEST_MODULE_INIT(SmithWatermanTest)

typedef BGIQD::ALIGN::SmithWaterman<char,std::string> Tester;
BGIQD::ALIGN::Schemes test_schemes{ 1 , 0 ,0 ,0 };
//
// by 
//  match = 1 
//  mismatch = 0 
//  insert = 0
//  delete = 0 
// 
// example 1
//      from
//          {"ABCBDAB"};
//          {"BDCABA"};
//      ==>
//          BCBA
//
// example 2
//      from
//          {"ACCGGTCGAGTGCGCGGAAGCCGGCCGAA"};
//          {"GTCGTTCGGAATGCCGTTGCTCTGTAAA"};
//      ==>
//          G T C G T C G G A A G C C G G
//


TEST(SWTest01)
{
    Tester test01 ;
    test01.schemes = test_schemes ;
    test01.ref = "ABCBDAB";
    test01.query = "BDCABA" ;
    test01.CheckLoadedData();
    test01.InitAfterDataLoaded();
    test01.FillMutrix();
    auto path = test01.GetResult();
    auto m1 = test01.AlignedElementsRef(path);
    auto m2 = test01.AlignedElementsRef(path);

    CHECK("BCBA" ,m1);
    CHECK("BCBA" ,m2);

}

TEST(SWTest02)
{
    Tester test01 ;
    test01.schemes = test_schemes ;
    test01.ref = "ACCGGTCGAGTGCGCGGAAGCCGGCCGAA";
    test01.query = "GTCGTTCGGAATGCCGTTGCTCTGTAAA" ;
    test01.CheckLoadedData();
    test01.InitAfterDataLoaded();
    test01.FillMutrix();
    auto path = test01.GetResult();
    auto m1 = test01.AlignedElementsRef(path);
    auto m2 = test01.AlignedElementsRef(path);

    CHECK("GTCGTCGGAAGCCGG" ,m1);
    CHECK("GTCGTCGGAAGCCGG" ,m2);
}
