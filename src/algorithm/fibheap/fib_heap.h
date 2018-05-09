#ifndef __ALGORITHM_FIBHEAP_H__
#define __ALGORITHM_FIBHEAP_H__

#include "algorithm/bi_list/bi_list.h"
#include <string>
#include <algorithm>
#include <vector>

namespace BGIQD {
    namespace FIBHEAP {

        template< class TKey , class TValue>
            struct Node
            {
                typedef TKey Key;

                typedef TValue Value ;

                typedef BGIQD::BILIST::BiList<Node> BiList;

                Node * father;

                Node * son;

                BiList list;

                Key    key;

                Value  value;

                bool   mark;

                int    degree ;

                void Init() 
                {
                    father = NULL ;
                    list.Init(this);
                    son = NULL ;
                    mark = false ;
                    degree = 0 ;
                }

                void AddChild(Node & child)
                {
                    assert( child.list.Single() ) ;
                    degree ++ ;
                    if( son == NULL )
                        son = &child ;
                    else
                    {
                        son->list.Insert(&child.list);
                        child.father = this ;
                    }
                }

                void RemoveChild( Node & child)
                {
                    assert( child.father == this );
                    degree -- ;
                    assert( degree >= 0 );
                    if( son == &child )
                    {
                        if( child.list.Single() ) 
                            son = NULL ;
                        else
                            son = child.list.left;
                    }
                    child.list.DeleteMe() ;
                }

                bool IsSingle() const 
                {
                    return list.Single() ;
                }

                Node & Next()
                {
                    return *(list.Forward()->self);
                }

                Node & Last()
                {
                    return *(list.Backword()->self);
                }

                void DeleteMe()
                {
                    list.DeleteMe();
                }

                void Insert( Node & node )
                {
                    list.Insert(&node.list);
                }
            };

        /**
         *
         * This FibHeap DO NOT care where nodes are stored !!!
         * It is user's duty to make sure all nodes are valid !!!
         * */
        template< class TKey
            , class TValue 
            , class TNode = Node< TKey , TValue>  >
            struct FibHeap
            {
                typedef TKey Key;
                typedef TValue Value ;
                typedef TNode Node;

                Node * min ;

                long long n ;

                FibHeap() : min( NULL ) , n (0) {}

                void Insert(Node & x)
                {
                    x.Init();

                    n ++ ;

                    if( min == NULL )
                    {
                        min = & x;
                        return ;
                    }

                    min->Insert(x);

                    if( x.key < min->key )
                    {
                        min = &x ;
                    }
                }

                Node & Min() const 
                {
                    if( min != NULL )
                        return *min ;
                    else
                        assert(0);
                }

                static FibHeap Union(const FibHeap & h1, const FibHeap &h2 )
                {
                    FibHeap H ;
                    H.min = h1.min ;
                    if( H.min != NULL && h2.min != NULL )
                        H.min->Insert(h2.min) ;
                    if( H.min == NULL || ( h2.min != NULL && h2.min->key < H.min->key ) )
                        H.min = h2.min ;
                    H.n = h1.n + h2.n ;
                }

                void DecreaseKey( Node & node , Key new_key)
                {
                    if ( node.key < new_key )
                    {
                        assert(0);
                        return ;
                    }

                    if( node.key == new_key )
                        return ;
                    node.key = new_key ;
                    if( node.father != NULL && node.key < node.father->key )
                    {
                        Cut(node,*node.father);
                        CascadingCut(*node.father);
                    }
                    assert( min );
                    if ( node.key < min->key )
                    {
                        min = & node ;
                    }
                }

                Node & ExtractMin()
                {
                    assert( min ) ;
                    auto & z = *min ;
                    if ( z.son != NULL )
                    {
                        Node * son = z.son ;
                        do
                        {
                            son->father = NULL ;
                            son = &(son->Next());
                        }
                        while( son != z.son );
                        min->Insert(*z.son);
                    }
                    if( z.IsSingle() )
                        min = NULL ;
                    else 
                    {
                        auto & next = z.Next() ;
                        z.DeleteMe() ;
                        min = &next ;
                        Consolidate();
                    }
                    n -- ;
                    return z;
                }

                bool Empty()
                {
                    return n == 0;
                }

                protected:

                    void Cut(Node & son , Node & father)
                    {
                        assert(min);
                        father.RemoveChild(son);
                        son.mark = false ;
                        son.father = NULL ;
                        min->Insert(son);
                    }

                    void CascadingCut(Node & father)
                    {
                        if( father.father == NULL )
                            return ;
                        auto & z = *father.father ;
                        if( father.mark == false )
                            father.mark = true ;
                        else 
                        {
                            Cut(father,z);
                            CascadingCut(z);
                        }
                    }
                    // make a big D(n) as 100
                    // 1.618^100 = 790,408,728,675,299,400,000
                    static const int N = 100 ;
                    struct DelayHeapLinkInfo
                    {
                        Node * y;
                        Node * x;
                    };

                    void Consolidate() 
                    {
                        assert( min ) ;
                        Node * A[N] ;
                        for( int i = 0 ; i < N ; i++ )
                            A[i] = NULL ;
                        if ( min->IsSingle() )
                            return ;
                        Node * a_son = min ;
                        Node * last = &(min->Last());

                        std::vector<DelayHeapLinkInfo> delayinfos;
                        do
                        {
                            Node * x = a_son ;
                            a_son = &(x->Next());
                            int d = x->degree ;
                            assert( d < N );
                            while( A[d] != NULL )
                            {
                                Node * y = A[d];
                                if ( x->key > y->key )
                                    std::swap(x,y);
                                //HeapLink(y,x);
                                DelayHeapLinkInfo tmp{ y , x } ;
                                delayinfos.emplace_back( tmp );
                                A[d] = NULL ;
                                d = d + 1 ;
                            }
                            assert( d < N );
                            A[d] = x ;
                        }
                        while( a_son != min );
                        for( size_t i = 0 ; i < delayinfos.size() ; i++ )
                        {
                            HeapLink(delayinfos[i].y , delayinfos[i].x);
                        }

                        for( int i = 0 ; i < N ; i ++ )
                        {
                            if( A[i] == NULL ) 
                                continue ;
                            if( A[i]->key < min->key )
                                min = A[i] ;
                        }

                    }

                    void HeapLink(Node * y , Node * x )
                    {
                        assert( x != NULL && y != NULL );
                        assert( x->father== NULL && y->father == NULL );
                        y->DeleteMe();
                        x->AddChild(*y);
                        y->mark = false ;
                    }
            };

    }

}

#endif //__ALGORITHM_FIBHEAP_H__
