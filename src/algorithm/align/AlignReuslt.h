#ifndef __ALGORITHM_ALIGN_ALIGNRESULT_H__
#define __ALGORITHM_ALIGN_ALIGNRESULT_H__

namespace BGIQD {
    namespace ALIGN {

        enum ResultType
        {
            Unknow = 0 ,
            LeftTop = 1 ,
            Left = 2 ,
            Top = 3 ,
        };

        struct AlignResult
        {
            int row_id ;
            int column_id ;
            ResultType type ;
        };
    }
}
#endif
