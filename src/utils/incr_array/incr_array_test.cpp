#include "utils/incr_array/incr_array.h"
#include "utils/unittest/Test.h"

TEST_MODULE_INIT(IncrArray)

typedef BGIQD::INCRARRAY::IncrArray<int> IntArray;

IntArray & GetTestData()
{
    static IntArray test(10);

    if( test.empty() )
    {
        for( int i = 0 ; i < 88 ; i++ )
            test.push_back(i);
    }

    return test ;
}

TEST(IncrArrayAcess)
{
    auto & test = GetTestData() ;

    for( int i = 0 ; i < 88 ; i++ )
    {
        CHECK( i , test[i]);
    }
    CHECK(88 ,test.size() );
    CHECK(90 ,test.capacity() );
    int j = 0 ;
    for( auto i = test.begin() ; i != test.end() ; i++ )
    {
        CHECK( j , *i);
        j++ ;
    }
    CHECK( 88 , j );
    j = 0 ;

    for( auto i = test.begin() ; i < test.end() ; i+= 2 )
    {
        CHECK( j , *i);
        j+= 2 ;
    }
    CHECK( 88 , j );

    j = 0; 
    for( auto & i : test )
    {
        CHECK( j , i);
        j ++ ;
    }
    CHECK( 88 , j );
}

TEST(IncrArrayClean)
{
    auto & test = GetTestData() ;
    test.clear();
    CHECK(0,test.size());
    CHECK(90,test.capacity());
}

TEST(IncrArrayDeepClean)
{
    auto & test = GetTestData() ;
    test.deep_clean();
    CHECK(0,test.size());
    CHECK(0,test.capacity());
}

TEST(IncrArrayResize)
{
    auto & test = GetTestData() ;
    test.resize(100);
    CHECK(100,test.size());
    test.resize(10);
    CHECK(10,test.size());
    CHECK(100,test.capacity());
}
