#ifndef __COMMON_MULTITHREAD_MAPREDUCE_H__
#define __COMMON_MULTITHREAD_MAPREDUCE_H__
#include "common/multithread/MultiThread.h"
#include <map>
#include <functional>
#include <vector>

namespace BGIQD{
namespace MultiThread{

template<class K2, class V2, class V3>
class Reducer
{
    public:
        typedef std::function< std::pair<K2 ,V3>(const K2 & ,const std::vector<V2>) > reducer;
        reducer func;
};

template<class K2, class V2>
class Combiner
{
    public:
        typedef std::function< std::pair<K2, std::vector<V2>> (const K2 & ,const std::vector<V2>) > combiner;
        combiner func;
};

// Not multi-thread safe. 
// A Shuffler must run in 1 thread , can be the same thread of mapper.
template<class K2, class V2>
class Shuffler
{
    public:
        typedef std::map<V2 ,std::vector<K2> > ShufflerResult;
        ShufflerResult && GetResults()
        {
            return std::move(m_results);
        }
        void AddResult(const K2& k ,const V2 & v )
        {
            m_results[k].emplace_back(v);
        }

    private:
        ShufflerResult m_results;
};

template<class K1 , class V1 , class K2 , class V2>
class Maper
{
    public:
        typedef std::function<std::pair<K2,V2>(K1 k , V1 v)> maper;
        maper func;
};


template<class K1 , class K2 , class V1 , class V2, class V3>
class MapReduce
{
    public:
        void SetMaper(int num , typename Maper<K1,V1,K2,V2>::maper func);
        void SetReducer(int num , typename Reducer<K2,V2,V3>::reducer func);
        void SetCombiner(int num, typename Combiner<K2,V2>::combiner func);
    private:
        std::vector<Maper<K1,V1,K2,V2>> m_mapper;
        std::vector<Shuffler<K2,V2>> m_shuffler;
        std::vector<Reducer<K2,V2,V3>> m_reducer;
        std::vector<Combiner<K2,V2>> m_combiner;

    private:
        MultiThread m_mapper_threads; //all result saved in shuffle
        MultiThread m_reducer_threads;
        MultiThread m_combiner_threads;

    public:
        

};

}
}

#endif //__COMMON_MULTITHREAD_MAPREDUCE_H__
