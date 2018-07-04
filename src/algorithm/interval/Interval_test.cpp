#include "algorithm/interval/Interval.h"
#include "common/test/Test.h"

TEST_MODULE_INIT(Interval)

typedef BGIQD::INTERVAL::Interval<int
    ,BGIQD::INTERVAL::IntervalType::Left_Close_Right_Close> IntervalInt01;

typedef BGIQD::INTERVAL::Interval<int
    ,BGIQD::INTERVAL::IntervalType::Left_Close_Right_Open> IntervalInt02;

typedef BGIQD::INTERVAL::Interval<int
    ,BGIQD::INTERVAL::IntervalType::Left_Open_Right_Close> IntervalInt03;

typedef BGIQD::INTERVAL::Interval<int
    ,BGIQD::INTERVAL::IntervalType::Left_Open_Right_Open> IntervalInt04;

TEST(TIntervalInt01)
{
    IntervalInt01 t(1,5);
    CHECK(true , t.IsContain(1));
    CHECK(true , t.IsContain(5));
    CHECK(true , t.IsContain(3));
    CHECK(false, t.IsContain(0));
    CHECK(false, t.IsContain(6));
}
TEST(TIntervalInt02)
{
    IntervalInt02 t(1,5);
    CHECK(true , t.IsContain(1));
    CHECK(false, t.IsContain(5));
    CHECK(true , t.IsContain(3));
    CHECK(false, t.IsContain(0));
    CHECK(false, t.IsContain(6));
}

TEST(TIntervalInt03)
{
    IntervalInt03 t(1,5);
    CHECK(false, t.IsContain(1));
    CHECK(true , t.IsContain(5));
    CHECK(true , t.IsContain(3));
    CHECK(false, t.IsContain(0));
    CHECK(false, t.IsContain(6));
}


TEST(TIntervalInt04)
{
    IntervalInt04 t(1,5);
    CHECK(false, t.IsContain(1));
    CHECK(false, t.IsContain(5));
    CHECK(true , t.IsContain(3));
    CHECK(false, t.IsContain(0));
    CHECK(false, t.IsContain(6));
}
