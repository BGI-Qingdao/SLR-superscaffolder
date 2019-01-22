#ifndef __ALGORITHM_ALIGN_SMITHWATERMAN_H__
#define __ALGORITHM_ALIGN_SMITHWATERMAN_H__

#include "algorithm/align/ScoreSchemes.h"
#include "algorithm/align/ScoreMatrix.h"
#include "algorithm/align/AlignReuslt.h"

#include <vector>
#include <stack>

namespace BGIQD {
    namespace ALIGN {

        template<class T , class V = std::vector<T> >
            struct SmithWaterman 
            {
                typedef T Item;
                typedef V DataContainer ;
                DataContainer ref;              // column
                DataContainer query;            // row
                BGIQD::ALIGN::Matrix the_matrix ;
                BGIQD::ALIGN::Schemes schemes;

                int  max_i ;
                int  max_j ;
                int max_value ;

                void CheckLoadedData()
                {
                    assert( schemes.match_score > 0 );
                    assert( schemes.insert_score <= 0 );
                    assert( schemes.delete_score<= 0 );
                    assert( schemes.mismatch_score<= 0 );
                    assert( ref.size() > 0 );
                    assert( query.size() > 0 );
                    return ;
                }

                void InitAfterDataLoaded()
                {
                    the_matrix.Init( ref.size() , query.size() );
                    for( int i = 0 ; i < the_matrix.Row() ; i++ )
                        the_matrix.Set(i,0,0);
                    for( int i = 0 ; i < the_matrix.Column() ; i++ )
                        the_matrix.Set(0,i,0);
                    max_i = 0 ;
                    max_j = 0 ;
                    max_value = 0 ;
                }

                void FillMutrix()
                {
                    for( int i = 1 ; i < the_matrix.Row() ; i++ )
                        for( int j = 1 ; j< the_matrix.Column() ; j++ )
                        {
                            int left_top =  the_matrix.Get(i-1,j-1);
                            int left  =  the_matrix.Get(i,j-1);
                            int top  =  the_matrix.Get(i-1,j);
                            Item & cs = ref.at(i) ;
                            Item & rs = query.at(i) ;
                            int from_lt = 0 ;
                            int from_l  = 0 ;
                            int from_t  = 0 ;
                            if ( cs == rs )
                                from_lt = left_top + schemes.match_score ;
                            else 
                                from_lt = left_top + schemes.mismatch_score ;
                            from_l = left + schemes.delete_score;
                            from_t = top + schemes.insert_score ;
                            int tmp = std::max( from_l , from_t );
                            tmp = std::max( tmp , from_lt ) ;
                            tmp = std::max( tmp , 0 );
                            the_matrix.Set( i , j , tmp );
                            if( tmp > max_value )
                            {
                                max_i = i ;
                                max_j = j ;
                            }
                        }
                }

                bool IsMatch( int threshold_score = 1 )
                {
                    return max_value >= threshold_score ;
                }

                DataContainer AlignedElementsRef( 
                        std::stack<BGIQD::ALIGN::AlignResult> result )
                {
                    DataContainer ret ;
                    while( ! result.empty() )
                    {
                        auto i = result.top() ;
                        if( i.type == BGIQD::ALIGN::ResultType::LeftTop )
                        {
                            ret.push_back( ref.at(i.column_id-1));
                        }
                        result.pop() ;
                    }
                    return ret ;
                }

                DataContainer AlignedElementsQuery( 
                        std::stack<BGIQD::ALIGN::AlignResult> result )
                {
                    DataContainer ret ;
                    while( ! result.empty() )
                    {
                        auto i = result.top() ;
                        if( i.type == BGIQD::ALIGN::ResultType::LeftTop )
                        {
                            ret.push_back( query.at(i.row_id-1));
                        }
                        result.pop() ;
                    }
                    return ret ;
                }

                std::stack<BGIQD::ALIGN::AlignResult> GetResult() const
                {
                    std::stack<BGIQD::ALIGN::AlignResult> ret ;
                    int curr = max_value ;
                    int i = max_i ;
                    int j = max_j ;
                    ret.push(
                            {max_i , max_j 
                            , BGIQD::ALIGN::ResultType::LeftTop} 
                            );
                    while( curr > 0 )
                    {
                        int left_top =  the_matrix.Get(i-1,j-1);
                        int left  =  the_matrix.Get(i,j-1);
                        int top  =  the_matrix.Get(i-1,j);
                        if( left >= top && left >= left_top )
                        {
                            j = j -1 ;
                            ret.push(
                                    {i , j 
                                    , BGIQD::ALIGN::ResultType::Left}
                                    );
                        }
                        else if ( top >= left && top >= left_top )
                        {
                            i = i -1 ;
                            ret.push(
                                    {i , j 
                                    , BGIQD::ALIGN::ResultType::Top}
                                    );
                        }
                        else
                        {
                            j = j -1 ;
                            i = i - 1 ;
                            ret.push(
                                    {i , j 
                                    , BGIQD::ALIGN::ResultType::LeftTop}
                                    );

                        }
                        curr = the_matrix.Get(i,j);
                    }
                    return ret ;
                }
            };
    }
}

#endif
