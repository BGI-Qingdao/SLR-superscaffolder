#ifndef __ALGORITHM_INTERVAL_INTERVAL_H__
#define __ALGORITHM_INTERVAL_INTERVAL_H__

#include <cassert>

namespace BGIQD{
    namespace INTERVAL{

        //struct IntervalType_Left_Close_Right_Close
        //{
        //    static bool IsLeftOpen() { return false ; }
        //    static bool IsRightOpen() { return false ; }
        //};


        enum IntervalType
        {
            Unknow = 0,
            Left_Close_Right_Close = 1 ,
            Left_Open_Right_Open = 2 ,
            Left_Close_Right_Open = 3 ,
            Left_Open_Right_Close = 4 ,
            Invalid = 5
        };

        template<class T , IntervalType t = IntervalType::Left_Close_Right_Close>
            struct Interval
            {
                typedef T ValueType;
                typedef IntervalType Type ;

                Interval() 
                {
                    type = t ;
                }

                Interval(const ValueType & l , const ValueType & r )
                {
                    type = t ;
                    min = l ;
                    max = r ;
                }

                Interval( const Interval & other )
                {
                    type = other.type ;
                    min = other.min ;
                    max = other.max ;
                }

                Interval & operator = ( const Interval & other )
                {
                    if( &other != this )
                    {
                        type = other.type ;
                        min = other.min ;
                        max = other.max ;
                    }
                    return *this;
                }

                Type type;
                ValueType min;
                ValueType max;

                bool IsContain(const ValueType& x )
                {
                    if( type == Type::Left_Close_Right_Close )
                    {
                        return x >= min && x <= max ;
                    }
                    else if( type == Type::Left_Close_Right_Open)
                    {
                        return x >= min && x < max ;
                    }
                    else if( type == Type::Left_Open_Right_Close)
                    {
                        return x > min && x <= max ;
                    }
                    else if( type == Type::Left_Open_Right_Open )
                    {
                        return x > min && x < max ;
                    }
                    assert(0);
                    return false ;
                }
            };
    }
}

#endif
