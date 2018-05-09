#include "common/test/Test.h"
#include "algorithm/fibheap/fib_heap.h"

TEST_MODULE_INIT(FibHeap)

typedef BGIQD::FIBHEAP::FibHeap<int , int > testFibHeap;
typedef BGIQD::FIBHEAP::Node<int , int > testFibHeapNode;
TEST(FibCommon)
{
    testFibHeapNode array[10];
    for( int i = 0 ; i < 10 ; i ++ )
    {
        array[i].key = std::rand() % 100 ;
        array[i].value = array[i].key * 2  ;
    }

    testFibHeap test;
    for( int i = 0 ; i < 10 ; i ++ )
    {
        test.Insert(array[i]);
    }
    int prev = -1 ;
    while( ! test.Empty() )
    {

        auto & min = test.ExtractMin() ;
        CHECK(true , (prev <= min.key )) ;
        std::cout<<min.key<<'\t'<<min.value<<std::endl;
        prev= min.key ;
    }
}



