#ifndef __ALGORITHM_BILIST_BILIST_H__
#define __ALGORITHM_BILIST_BILIST_H__

#include <cassert>
#include <string>
namespace BGIQD{
    namespace BILIST{

        template<class T>
            struct BiList
            {
                typedef BiList * BiListPtr ;
                typedef T *      ContainerPtr;
                BiListPtr left;
                BiListPtr right;
                const T * self;

                void Init(const T * s) 
                {
                    left = this;
                    right = this;
                    self = s ;
                }

                void Insert( BiListPtr node)
                {
                    assert( node != NULL );
                    assert( left != NULL && right != NULL );
                    assert( node->left != NULL && node->right != NULL );

                    auto A = this ;
                    auto B = node ;
                    auto ar = A->right ;
                    auto br = B->right ;
                    A->right->left = B ;
                    B->right->left = A ;
                    A->right = br ;
                    B->right= ar ;

                }

                void DeleteMe() 
                {
                    assert( left != NULL && right != NULL );
                    if( left != this  && right != this )
                    {
                        left->right = right ;
                        right->left = left ;
                    }
                    else
                    {
                        assert( left == this && right == this );
                    }
                    Init(self);
                }
            };
    }
}
#endif //__ALGORITHM_BILIST_BILIST_H__
