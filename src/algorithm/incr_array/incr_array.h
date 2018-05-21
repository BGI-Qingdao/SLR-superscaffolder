#ifndef __ALGORITHM_INCR_ARRAY_INCR_ARRAY_H__
#define __ALGORITHM_INCR_ARRAY_INCR_ARRAY_H__

#include <cstddef>
#include <cassert>

#include <vector>
#include <iterator>
//
// Incr Array .
//
// A substitute of vector which has below feather
//
// good way
//  * DO NOT copy old element while expend size.
//  * element pointer ALWAYS valid.
// bad way
//  * new element can ONLY append at the end of array.
//  * DO NOT suppert earse element.
//  * DO NOT suppert copy-construct and assign-construct.
//      * can support but why need this ?
//

namespace BGIQD {
    namespace INCRARRAY {

        template<class T>
            struct iterator ;

        template<class T>
            struct  IncrArray
            {
                typedef T Element;

                typedef struct iterator<IncrArray> iterator;

                virtual ~IncrArray()
                {
                    for( auto  & i : m_headers )
                    {
                        if( i != NULL )
                        {
                            delete [] i ;
                            i = NULL ;
                        }
                    }
                }

                IncrArray(size_t block_size)
                {
                    m_block_size = block_size ;
                    m_curr = 0 ;
                    m_headers.push_back(new Element[m_block_size]);
                    m_capacity = m_block_size;
                }

                IncrArray( const IncrArray & );

                IncrArray & operator = ( const IncrArray & ) ;

                void push_back(const Element & e)
                {
                    assert( m_curr <= m_capacity );
                    if ( m_curr ==  m_capacity )
                    {
                        m_headers.push_back(new Element[m_block_size]);
                        m_capacity += m_block_size ;
                    }
                    assert( m_curr < m_capacity );
                    operator[](m_curr++) = e ;
                }

                Element & operator[](size_t i )
                {
                    assert( i < m_curr );
                    size_t block_num = i / m_block_size ;
                    size_t block_i= i % m_block_size ;
                    if( block_num >= m_headers.size() )
                    {
                        assert(0);
                    }
                    return m_headers[block_num][block_i];
                }

                iterator begin()
                {
                    iterator ret(*this,0);
                    return ret ;
                }

                iterator end()
                {
                    iterator ret(*this,m_curr);
                    return ret ;
                }

                bool empty() const { return m_curr == 0 ;}

                void clear() { m_curr = 0 ; }

                size_t size() const { return m_curr ; }

                size_t capacity() const { return m_curr ; }

                void resize( size_t i ) 
                {
                    if( i <= m_capacity ) 
                        return ;
                    do {
                        m_headers.push_back(new Element[m_block_size]);
                        m_capacity += m_block_size ;
                    }while( m_capacity < i );
                }

                protected:
                size_t m_block_size ;
                size_t m_curr ;
                size_t m_capacity;
                std::vector<Element * > m_headers;
            };

        template<class T>
            struct iterator : public std::iterator<typename T::Element
                              , size_t 
                              , typename T::Element * 
                              , typename T::Element &
                              , std::random_access_iterator_tag>
        {
            typedef T Base;

            typedef typename T::Element * pointer;

            typedef typename T::Element & reference;

            iterator( Base & b , size_t i ) : base(&b) , curr(i) {}

            iterator( const iterator & o ) : base(o.base) , curr(o.curr) {}

            iterator & operator = ( const iterator & o )
            {
                base = o.base ;
                curr = o.curr ;
                return *this ;
            }

            virtual ~iterator() {} ;

            bool operator <= ( const iterator & o )const { 
                assert ( base == o.base );
                return curr <= o.curr ; 
            }

            bool operator >= ( const iterator & o ) const { 
                assert ( base == o.base );
                return curr >= o.curr ; 
            }

            bool operator != ( const iterator & o ) const {
                assert ( base == o.base );
                return curr == o.curr ; 
            }
            bool operator < ( const iterator & o )const { 
                assert ( base == o.base );
                return curr < o.curr ; 
            }

            bool operator > ( const iterator & o ) const { 
                assert ( base == o.base );
                return curr > o.curr ; 
            }

            bool operator == ( const iterator & o ) const {
                assert ( base == o.base );
                return curr == o.curr ; 
            }

            reference operator*() 
            {
                return (*base)[curr];
            }

            pointer operator->()
            {
                return &((*base)[curr]);
            }

            iterator & operator++ (int i)
            {
                curr += i ;
                return *this;
            }

            iterator & operator-- (int i)
            {
                curr -= i ;
                return *this;
            }
            iterator & operator++ ()
            {
                curr ++ ;
                return *this;
            }

            iterator & operator-- ()
            {
                curr -- ;
                return *this;
            }

            iterator & operator += ( size_t i)
            {
                curr += i ;
                return *this ;
            }

            iterator & operator -= ( size_t i)
            {
                curr += i ;
                return *this ;
            }

            iterator  operator + ( size_t i) const 
            {
                iterator ret(*this);
                ret.curr += i ;
                return ret;
            }

            iterator  operator - ( size_t i) const 
            {
                iterator ret(*this);
                ret.curr -= i ;
                return ret;
            }

            protected:
            Base * base;
            size_t curr ;
        };

    } // INCRARRAY
} // BGIQD

#endif //__ALGORITHM_INCR_ARRAY_INCR_ARRAY_H__
