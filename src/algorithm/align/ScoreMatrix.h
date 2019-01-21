#ifndef __ALGORITHM_ALIGN_SCOREMATRIX_H__
#define __ALGORITHM_ALIGN_SCOREMATRIX_H__

namespace BGIQD {
    namespace ALIGN {

        struct  Matrix
        {
            void Init( int l1 , int l2 ) 
            {
                
            }
            int Len1() const  { return len1 +1 ;}
            int Len2() const  { return len2 +1 ;}

            int Get( int x , int y ) const 
            {
                return  
            }
            private:
            int len1;
            int len2;
            std::vector<int> data;
        };
    }
}

#endif
