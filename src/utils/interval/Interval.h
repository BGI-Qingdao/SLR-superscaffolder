#ifndef __ALGORITHM_INTERVAL_INTERVAL_H__
#define __ALGORITHM_INTERVAL_INTERVAL_H__

#include <cassert>
#include <algorithm>
#include <string>

/**********************************************************
 *
 * @Brief :
 *   Define struct for interval like [1,7] or (10,20)
 *   Implement interface to detect contain or overlap between 2 interval.
 *   Implement interface to detect contain between one point and one interval.
 *  
 * *******************************************************/
namespace BGIQD{
    namespace INTERVAL{

        // Define traits to handle [] [) (] and ()
        enum IntervalType
        {
            Unknow = 0,
            Left_Close_Right_Close = 1 ,
            Left_Open_Right_Open = 2 ,
            Left_Close_Right_Open = 3 ,
            Left_Open_Right_Close = 4 ,
            Invalid = 5
        };

        // extract traits
        template<IntervalType type>
            struct IntervalInfo
            {
            };

        template<>
            struct IntervalInfo<Left_Close_Right_Open>
            {
                static const  bool LeftOpen = false;
                static const bool RightOpen = true;
            };
        template<>
            struct IntervalInfo<Left_Open_Right_Close>
            {
                static const  bool LeftOpen = true;
                static const bool RightOpen = false;
            };

        template<>
            struct IntervalInfo<Left_Open_Right_Open>
            {
                static const  bool LeftOpen = true;
                static const bool RightOpen = true;
            };
        template<>
            struct IntervalInfo<Left_Close_Right_Close>
            {
                static const  bool LeftOpen = false;
                static const bool RightOpen = false ;
            };

        // The interval class
        template<class T , IntervalType t = IntervalType::Left_Close_Right_Close>
            struct Interval
            {
                typedef T ValueType;
                typedef IntervalType Type ;

                static const Type type = t;
                ValueType min;
                ValueType max;

                Interval() 
                {
                }

                Interval(const ValueType & l , const ValueType & r )
                {
                    min = l ;
                    max = r ;
                }

                Interval( const Interval & other )
                {
                    min = other.min ;
                    max = other.max ;
                }

                Interval & operator = ( const Interval & other )
                {
                    if( &other != this )
                    {
                        min = other.min ;
                        max = other.max ;
                    }
                    return *this;
                }
                // ONLY USED FOR CONTAINER !!!
                bool operator < ( const Interval & other ) const
                {
                    return min < other.min ;
                }
                bool operator == ( const Interval & other ) const
                {
                    return min == other.min && max == other.max ;
                }
                std::string ToString() const 
                {
                    std::string ret ;
                    if ( IntervalInfo<type>::LeftOpen )
                    {
                        ret += "( ";
                    }
                    else
                    {
                        ret += "[ ";
                    }
                    ret += std::to_string(min);
                    ret += " , ";
                    ret += std::to_string(max);
                    if ( IntervalInfo<type>::RightOpen )
                    {
                        ret += " )";
                    }
                    else
                    {
                        ret += " ]";
                    }

                    return ret;
                }

                // Is a point contained by this interval ?
                bool IsContain(const ValueType& x ) const
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

                ValueType Len() const 
                {
                    return max - min ;
                }

                // Is an other interval contained by this interval ?
                bool IsContain( const Interval & other ) const 
                {
                    if( min > other.min )
                        return false ;
                    if( max < other.max )
                        return false ;
                    return true ;
                }

                Interval Overlap( const Interval & other )
                {
                    Interval ret ;
                    ret.max = 0 ;
                    ret.min = 0 ;
                    // Check contain
                    if( IsContain(other) )
                        return other;
                    if( other.IsContain( *this) )
                        return *this ;
                    // Check no overlap
                    if( ( !IntervalInfo<type>::RightOpen)  && (! IntervalInfo<type>::LeftOpen ))
                    {
                        if( max < other.min || min > other.max )
                            return ret ;
                    }
                    else
                    {
                        if( max <= other.min || min >= other.max) 
                            return ret ;
                    }
                    // Return overlap
                    ret.min = std::max(min,other.min);
                    ret.max = std::min(max,other.max);
                    return ret ;
                }
            };
    }
}

#endif
